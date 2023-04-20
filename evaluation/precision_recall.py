#precision/recall vs GIAB
import hail as hl
import pyspark
from wes_qc.utils.utils import parse_config


def prepare_alspac_htfile(mtfile: str, rf_htfile: str, mtdir: str) -> hl.Table:
    '''
    Filter mtfile to GIAB sample only. Annotate with rf bin, genotype hard filters, remove unused variants.
    :param str mtfile: Matrixtable file after splitting
    :param str rf_htfile: Random forest ht file
    :param str mtdir: MatrixTable directory
    :return: hl.Table
    '''
    mt = hl.read_matrix_table(mtfile)
    #filter to GIAB sample of interest
    sample = 'EGAN00003332049'#GIAB12878/HG001
    mt = mt.filter_cols(mt.s == sample)

    #hard filters
    filter_condition = (
        (mt.GT.is_het() & (mt.HetAB < 0.3)) | 
        (mt.DP < 10) |
        (mt.GQ < 10)
    )
    mt = mt.annotate_entries(
        hard_filters = hl.if_else(filter_condition, 'Fail', 'Pass')
    )
    mt = mt.filter_entries(mt.hard_filters == 'Pass')

    #filter to autosomes
    mt = mt.filter_rows(mt.locus.in_autosome())

    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(mt.variant_qc.n_non_ref > 0)

    rf_ht = hl.read_table(rf_htfile)
    # keep only vars with rank_id = rank
    rf_ht = rf_ht.filter(rf_ht.rank_id == 'rank')
    # annotate mt with score and bin
    mt = mt.annotate_rows(
        info=mt.info.annotate(
            rf_score=rf_ht[mt.row_key].score)
    )
    mt = mt.annotate_rows(
        info=mt.info.annotate(
            rf_bin=rf_ht[mt.row_key].bin)
    )

    alspac_vars = mt.rows()
    tmpht = mtdir + "tmp.ht"
    alspac_vars = alspac_vars.checkpoint(tmpht, overwrite = True)
    
    return alspac_vars


def prepare_giab_ht(giab_vcf: str, mtdir: str) -> hl.Table:
    '''
    Get GIAB ht from vcf file
    :param str giab_vcf: path of input VCF file
    :param str mtdir: MatrixTable directory
    :return: hl.Table
    '''
    mt = hl.import_vcf(giab_vcf, force_bgz = True, reference_genome='GRCh38')
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(mt.variant_qc.n_non_ref > 0)

    giab_vars = mt.rows()

    tmpht = mtdir + "tmpg.ht"
    giab_vars = giab_vars.checkpoint(tmpht, overwrite = True)

    return giab_vars


def get_precision_recall(giab_vars: hl.Table, alspac_vars: hl.Table, mtdir: str) -> tuple:
    '''
    Get precision and recall for two sets of variants, reference set (GIAB) and test set (ALSPAC)
    :param hl.Table giab_vars: GIAB variants (reference set)
    :param hl.Table alspac_vars: ALSPAC variants (test set)
    :param str mtdir: MatrixTable directory
    :return: tuple
    '''
    #checkpoint datasets to temporary files for speed
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


def calculate_precision_recall(alspac_ht: hl.Table, giab_ht: hl.Table, mtdir: str) -> dict:
    '''
    Calculate precision and recall for each RF bin for both SNPs and indels
    :param hl.Table alspace_ht: ALSPAC variants hail Table
    :param hl.Table giab_ht: GIAB variants hail Table
    :return: dict
    :param str mtdir: MatrixTable directory
    '''
    results = {'snv':{}, 'indel':{}}
    n_bins = 102 #there are 101 bins but use 102 for range to work

    giab_snvs = giab_ht.filter(hl.is_snp(giab_ht.alleles[0], giab_ht.alleles[1]))
    giab_indels = giab_ht.filter(hl.is_indel(giab_ht.alleles[0], giab_ht.alleles[1]))
    giab_indels = giab_indels.key_by(**hl.min_rep(giab_indels.locus, giab_indels.alleles))

    for bin in range(1,n_bins):
        results['snv'][bin] = {}
        results['indel'][bin] = {}
        alspac_filtered_ht = alspac_ht.filter(alspac_ht.info.rf_bin <= bin)
        alspac_var_count = alspac_filtered_ht.count()
        print("bin = " + str(bin) + " n vars = " + str(alspac_var_count))#for sanity check
        alspac_snvs = alspac_filtered_ht.filter(hl.is_snp(alspac_filtered_ht.alleles[0], alspac_filtered_ht.alleles[1]))
        alspac_indels = alspac_filtered_ht.filter(hl.is_indel(alspac_filtered_ht.alleles[0], alspac_filtered_ht.alleles[1]))
        alspac_indels = alspac_indels.key_by(**hl.min_rep(alspac_indels.locus, alspac_indels.alleles))

        snv_prec, snv_recall, snv_tp, snv_fp, snv_fn = get_precision_recall(giab_snvs, alspac_snvs, mtdir)
        indel_prec, indel_recall, indel_tp, indel_fp, indel_fn = get_precision_recall(giab_indels, alspac_indels, mtdir)

        results['snv'][bin]['precision'] = snv_prec
        results['snv'][bin]['recall'] = snv_recall
        results['indel'][bin]['precision'] = indel_prec
        results['indel'][bin]['recall'] = indel_recall
        results['snv'][bin]['TP'] = snv_tp
        results['snv'][bin]['FP'] = snv_fp
        results['snv'][bin]['FN'] = snv_fn
        results['indel'][bin]['TP'] = indel_tp
        results['indel'][bin]['FP'] = indel_fp
        results['indel'][bin]['FN'] = indel_fn

    return results


def print_results(results: dict, plot_dir: str):
    '''
    print results to file
    :param dict results: results dict
    :param str plotdir: Output plot directory
    '''
    outfile = plot_dir + "/precision_recall_dp10_gq10_ab03.txt"
    n_bins = 102
    header = ("\t").join(['bin', 'snv_precision', 'snv_recall', 'indel_precision', 'indel_recall', 'snv_tp', 'snv_fp',
        'snv_fn', 'indel_tp', 'indel_fp', 'indel_fn'])
    with open(outfile, 'w') as o:
        o.write(header)
        o.write("\n")
        for bin in range(1,n_bins):
            binstr = str(bin)
            snv_prec = str(results['snv'][bin]['precision'])
            snv_recall = str(results['snv'][bin]['recall'])
            indel_prec = str(results['indel'][bin]['precision'])
            indel_recall = str(results['indel'][bin]['recall'])
            snv_tp = str(results['snv'][bin]['TP'])
            snv_fp = str(results['snv'][bin]['FP'])
            snv_fn = str(results['snv'][bin]['FN'])
            indel_tp = str(results['indel'][bin]['TP'])
            indel_fp = str(results['indel'][bin]['FP'])
            indel_fn = str(results['indel'][bin]['FN'])
            outline = ("\t").join([binstr, snv_prec, snv_recall, indel_prec, indel_recall, snv_tp, snv_fp, snv_fn,
                indel_tp, indel_fp, indel_fn])
            o.write(outline)
            o.write("\n")
        

def main():
    # set up
    inputs = parse_config()
    mtdir = inputs['matrixtables_lustre_dir']
    rf_dir = inputs['var_qc_rf_dir']
    resourcedir = inputs['resource_dir']
    plot_dir = inputs['plots_dir_local']

    # initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")


    rf_htfile = rf_dir + "fef0f446" + "/_gnomad_score_binning_tmp.ht"
    mtfile = mtdir + "mt_varqc_splitmulti.mt"

    alspac_vars_ht = prepare_alspac_htfile(mtfile, rf_htfile, mtdir)

    giab_vcf = resourcedir + "HG001_GRCh38_benchmark.interval.illumina.vcf.gz"

    giab_ht = prepare_giab_ht(giab_vcf, mtdir)

    results = calculate_precision_recall(alspac_vars_ht, giab_ht, mtdir)

    print_results(results, plot_dir)


if __name__ == '__main__':
    main()