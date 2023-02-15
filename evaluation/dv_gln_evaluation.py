#evaluate dv/gln filtered calls
import hail as hl
import pyspark
import datetime
from wes_qc.utils.utils import parse_config
from wes_qc.evaluation.variant_counts_per_cq_controls.annotation_functions import annotate_cq, annotate_gnomad, get_counts_per_cq, get_ca_fractions

def get_variant_counts_per_cq_and_t_u(mtfile, cqfile, mtdir, pedfile, gnomad_htfile) -> float:
    '''
    :param str mtfile: ALSPAC mt file
    :param str cqfile: Most severe consequence annotation from VEP
    :param str mtdir: matrixtable directory
    :param str pedfile: pedigree file
    :param str gnomad_htfile: gnomad hail table file
    :return: flaot
    '''
    mt = hl.read_matrix_table(mtfile)
    #annotate mt with consequences and gnomad
    mt_cq = annotate_cq(mt, cqfile)
    mt_cq_gnomad = annotate_gnomad(mt_cq, gnomad_htfile)
    cq_outfile = mtdir + "consequences_counts.txt"
    ca_outfile = mtdir + "ca_counts.txt"
    get_counts_per_cq(mt_cq_gnomad, cq_outfile)
    get_ca_fractions(mt_cq_gnomad, ca_outfile)

    pedigree = hl.Pedigree.read(pedfile)
    #list of samples in trios
    trio_sample_ht = hl.import_fam(pedfile)
    sample_list = trio_sample_ht.id.collect() + trio_sample_ht.pat_id.collect() + trio_sample_ht.mat_id.collect()

    tu_ratio = get_trans_untrans(mt_cq_gnomad, pedigree, sample_list, mtdir)

    return tu_ratio
    

def get_trans_untrans(mt: hl.MatrixTable, pedigree: hl.Pedigree, sample_list: list, mtdir: str) -> float:
    '''
    get transmitted/untransmitted ratio
    :param hl.MatrixTable mt: matrixtable
    :param hl.Pedigree pedigree: Hail Pedigree
    :param list sample_list: List of samples in trios
    :param str mtdir: matrixtable directory
    :return float:
    '''
    #filter to synonymous
    mt_syn = mt.filter_rows(mt.consequence == 'synonymous_variant')
        #restrict to samples in trios, annotate with AC and filter to trio AC == 1 or 2
    mt2 = mt_syn.filter_cols(hl.set(sample_list).contains(mt_syn.s))
    mt2 = hl.variant_qc(mt2, name='varqc_trios')
    tmpmt3 = mtdir + "tmp3x.mt"
    mt2 = mt2.checkpoint(tmpmt3, overwrite = True)
    #split to potentially transitted/untransmitted
    untrans_mt = mt2.filter_rows(mt2.varqc_trios.AC[1] == 1)
    tmpmt4 = mtdir + "tmp4x.mt"
    untrans_mt = untrans_mt.checkpoint(tmpmt4, overwrite = True)
    trans_mt = mt2.filter_rows(mt2.varqc_trios.AC[1] == 2)
    tmpmt5 = mtdir + "tmp5x.mt"
    trans_mt = trans_mt.checkpoint(tmpmt5, overwrite = True)
        #run tdt function for potential trans and untrans
    tdt_ht_trans = hl.transmission_disequilibrium_test(trans_mt, pedigree)
    tdt_ht_untrans = hl.transmission_disequilibrium_test(untrans_mt, pedigree)
    trans_sing = tdt_ht_trans.filter((tdt_ht_trans.t == 1) & (tdt_ht_trans.u == 0))
    trans = trans_sing.count()
    untrans_sing = tdt_ht_untrans.filter((tdt_ht_untrans.t == 0) & (tdt_ht_untrans.u == 1))
    untrans = untrans_sing.count()
    if untrans > 0:
        ratio = trans/untrans
    else:
        ratio = 0

    return ratio


def precision_recall(mtfile: str, giab_vcffile: str, mtdir: str) -> dict:
    '''
    :param str mtfile: ALSPAC mt file
    :param str giab_vcffile: GIAB VCF file
    :param str mtdir: matrixtable directory
    :return: dict containing precision and recall values
    '''

    mt = hl.read_matrix_table(mtfile)
    #filter to GIAB sample and create HT of vars in this sample
    sample = 'EGAN00003332049'#GIAB12878/HG001
    mt = mt.filter_cols(mt.s == sample)
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(mt.variant_qc.n_non_ref > 0)
    alspac_vars = mt.rows()
    tmpht = mtdir + "tmp.ht"
    alspac_vars = alspac_vars.checkpoint(tmpht, overwrite = True)
    alspac_snvs = alspac_vars.filter(hl.is_snp(alspac_vars.alleles[0], alspac_vars.alleles[1]))
    alspac_indels = alspac_vars.filter(hl.is_indel(alspac_vars.alleles[0], alspac_vars.alleles[1]))
    alspac_indels = alspac_indels.key_by(**hl.min_rep(alspac_indels.locus, alspac_indels.alleles))

    #create giab HT
    giabmt = hl.import_vcf(giab_vcffile, force_bgz = True, reference_genome='GRCh38')
    giabmt = hl.variant_qc(giabmt)
    giabmt = giabmt.filter_rows(giabmt.variant_qc.n_non_ref > 0)
    giab_vars = giabmt.rows()
    tmpht = mtdir + "tmpg.ht"
    giab_vars = giab_vars.checkpoint(tmpht, overwrite = True)
    giab_snvs = giab_vars.filter(hl.is_snp(giab_vars.alleles[0], giab_vars.alleles[1]))
    giab_indels = giab_vars.filter(hl.is_indel(giab_vars.alleles[0], giab_vars.alleles[1]))
    giab_indels = giab_indels.key_by(**hl.min_rep(giab_indels.locus, giab_indels.alleles))

    #get precision recall for SNVs and indels
    snv_prec, snv_recall, snv_tp, snv_fp, snv_fn = get_precision_recall(giab_snvs, alspac_snvs, mtdir)
    indel_prec, indel_recall, indel_tp, indel_fp, indel_fn = get_precision_recall(giab_indels, alspac_indels, mtdir)

    results = {
        'snv_prec': snv_prec,
        'snv_recall': snv_recall,
        'indel_prec': indel_prec,
        'indel_recall': indel_recall
    }

    return results


def get_precision_recall(giab_vars: hl.Table, alspac_vars: hl.Table, mtdir: str) -> tuple:
    '''
    Get precision and recall for two sets of variants, reference set (GIAB) and test set (ALSPAC)
    :param hl.Table giab_vars: GIAB variants (reference set)
    :param hl.Table alspac_vars: ALSPAC variants (test set)
    :param str mtdir: MatrixTable directory
    :return: tuple
    '''

    tmpht1 = mtdir + "tmp1.ht"
    tmpht2 = mtdir + "tmp2.ht"
    giab_vars = giab_vars.checkpoint(tmpht1, overwrite = True)
    alspac_vars = alspac_vars.checkpoint(tmpht2, overwrite = True)

    print("get intersects")
    vars_in_both = giab_vars.semi_join(alspac_vars)
    giab_only = giab_vars.anti_join(alspac_vars)
    alspac_only = alspac_vars.anti_join(giab_vars)
    print("count_vars")
    tp = vars_in_both.count()
    fn = giab_only.count()
    fp = alspac_only.count()

    precision = tp / (tp + fp)
    recall = tp / (tp + fn)

    return precision, recall, tp, fp, fn


def main():
    #set up
    inputs = parse_config()
    mtdir = inputs['matrixtables_lustre_dir']
    resourcedir = inputs['resource_dir']

    # initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    mtfile_filtered = mtdir + "mt_after_hard_filter_gt.mt"
    giab_vcf = resourcedir + "HG001_GRCh38_benchmark.interval.illumina.vcf.gz"
    cqfile =  "file:///lustre/scratch123/qc/glnexus_vep/all_consequences_dv_gln.txt"
    pedfile = resourcedir + "trios.ped"
    gnomad_htfile = resourcedir + "gnomad.exomes.r2.1.1.sites.liftover_grch38.ht"

    tu_ratio = get_variant_counts_per_cq_and_t_u(mtfile_filtered, cqfile, mtdir, pedfile, gnomad_htfile)
    precison_recall = get_precision_recall(mtfile_filtered, giab_vcf, mtdir)

    print("transmitted/untransmitted ration for synonymous singletons: " + str(tu_ratio))
    print("Precision/recall:")
    print("Precision SNP: " + str(precision_recall['snv_prec']))
    print("Recall SNP: " + str(precision_recall['snv_recall']))
    print("Precision indel: " + str(precision_recall['indel_prec']))
    print("Recall indel: " + str(precision_recall['snv_recall']))


if __name__ == '__main__':
    main()