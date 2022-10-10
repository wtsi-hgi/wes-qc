#precision/recall vs GIAB
import hail as hl
import pyspark
from wes_qc.utils.utils import parse_config


def prepare_alspac_htfile(mtfile: str, rf_htfile: str) -> hl.Table:
    '''
    Filter mtfile to GIAB sample only. Annotate with rf bin, genotype hard filters, remove unused variants.
    :param str mtfile: Matrixtable file after splitting
    :param str rf_htfile: Random forest ht file
    :return: hl.Table
    '''
    mt = hl.read_matrix_table(mtfile)
    #filter to GIAB sample of interest
    sample = 'EGAN00003332049'#GIAB12878/HG0001
    mt = mt.filter_cols(mt.s == sample)

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
    
    return alspac_vars


def prepare_giab_ht(giab_vcf: str) -> hl.Table:
    '''
    Get GIAB ht from vcf file
    :param str giab_vcf: path of input VCF file
    :return: hl.Table
    '''
    mt = hl.import_vcf(giab_vcf, force_bgz = True, reference_genome='GRCh38')
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(mt.variant_qc.n_non_ref > 0)

    giab_vars = mt.rows()
    return giab_vars


def get_precision_recall(giab_vars, alspac_vars) -> tuple:
    '''
    Get precision and recall for two sets of variants, reference set (GIAB) and test set (ALSPAC)
    :param hl.Table giab_vars: GIAB variants (reference set)
    :param hl.Table alspac_vars: ALSPAC variants (test set)
    :return: tuple
    '''
    vars_in_both = giab_vars.semi_join(alspac_vars)
    giab_only = giab_vars.anti_join(alspac_vars)
    alspac_only = alspac_vars.anti_join(giab_vars)

    TP = vars_in_both.count()
    FN = giab_only.count()
    FP = alspac_only.count()

    precision = TP/(TP + FP)
    recall = TP / (TP + FN)

    return precision, recall


def calculate_precision_recall(alspac_ht, giab_ht) -> dict:
    '''
    Calculate precision and recall for each RF bin for both SNPs and indels
    :param hl.Table alspace_ht: ALSPAC variants hail Table
    :param hl.Table giab_ht: GIAB variants hail Table
    :return: dict 
    '''
    results = {'snv':{}, 'indel':{}}
    n_bins = 102 #there are 101 bins but use 102 for range to work

    giab_snvs = giab_ht.filter(hl.is_snp(giab_ht.alleles[0], giab_ht.alleles[1]))
    giab_indels = giab_ht.filter(hl.is_indel(giab_ht.alleles[0], giab_ht.alleles[1]))

    for bin in range(1,n_bins):
        results['snv'][bin] = {}
        results['indel'][bin] = {}
        alspac_filtered_ht = alspac_ht.filter(alspac_ht.info.bin <= bin)
        alspac_var_count = alspac_filtered_ht.count()
        print("bin = " + str(bin) + " n vars = " + str(alspac_var_count))#for sanity check
        alspac_snvs = alspac_filtered_ht.filter(hl.is_snp(alspac_filtered_ht.alleles[0], alspac_filtered_ht.alleles[1]))
        alspac_indels = alspac_filtered_ht.filter(hl.is_indel(alspac_filtered_ht.alleles[0], alspac_filtered_ht.alleles[1]))

        snv_prec, snv_recall = get_precision_recall(giab_snvs, alspac_snvs)
        indel_prec, indel_recall = get_precision_recall(giab_indels, alspac_indels)

        results['snv'][bin]['precision'] = snv_prec
        results['snv'][bin]['recall'] = snv_recall
        results['indel'][bin]['precision'] = indel_prec
        results['indel'][bin]['recall'] = indel_recall

    return results


def main():
    # set up
    inputs = parse_config()
    mtdir = inputs['matrixtables_lustre_dir']
    rf_dir = inputs['var_qc_rf_dir']
    resourcedir = inputs['resource_dir']

    # initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")


    rf_htfile = rf_dir + "6617f838" + "/_gnomad_score_binning_tmp.ht"
    mtfile = mtdir + "mt_varqc_splitmulti.mt"

    alspac_vars_ht = prepare_alspac_htfile(mtfile, rf_htfile)

    giab_vcf = resourcedir + "HG001_GRCh38_benchmark.interval.illumina.vcf.gz"

    giab_ht = prepare_giab_ht(giab_vcf)

    results = calculate_precision_recall(alspac_vars_ht, giab_ht)
    print(results)


if __name__ == '__main__':
    main()