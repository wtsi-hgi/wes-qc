#get variants in GIAB mthat fail ALSPAC genotype filters
import hail as hl
import pyspark
from wes_qc.utils.utils import parse_config


def annotate_mtfile_rf(mtfile: str, rf_htfile: str) -> hl.MatrixTable:
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

    return mt


def apply_hard_filters(mt: hl.MatrixTable) -> hl.MatrixTable:
    '''
    Apply genotype hard filters
    '''
    #hard filters
    filter_condition = (
        (mt.GT.is_het() & (mt.HetAB < 0.3)) | 
        (mt.DP < 10) |
        (mt.GQ < 20)
    )
    mt = mt.annotate_entries(
        hard_filters = hl.if_else(filter_condition, 'Fail', 'Pass')
    )

    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(mt.variant_qc.n_non_ref > 0)

    return mt


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

    tmpht = mtdir + "tmpg1.ht"
    giab_vars = giab_vars.checkpoint(tmpht, overwrite = True)

    return giab_vars


def get_vars_in_giab_fail_hard_filters(alspac_mt_giab_sample_hard_filters: hl.MatrixTable, giab_ht: hl.Table, mtout: str, mtdir: str):
    '''
    Get ALSPAC variants which are in GIAB but fail hard filters in ALSPAC
    :param hl.MatrixTable alspac_mt_giab_sample_hard_filters: Matrixtable of ALSPAC data for GIAB sample with hard filter annotation
    :param hl.Table giab_ht: Hail Table of GIAB variants
    :param str mtout: path to output matrixtable file
    :param str mtdir: path to matrixtable directory
    '''
    tmpmt = mtdir + "tmp1.mt"
    alspac_mt_giab_sample_hard_filters = alspac_mt_giab_sample_hard_filters.checkpoint(tmpmt, overwrite = True)
    mt = alspac_mt_giab_sample_hard_filters.semi_join_rows(giab_ht)

    tmpmt2 = mtdir + "tmp2.mt"
    mt = mt.checkpoint(tmpmt2, overwrite = True)

    mtfails = mt.filter_entries(mt.hard_filters == 'Fail')

    mtfails = hl.variant_qc(mtfails)
    mtfails = mtfails.filter_rows(mtfails.variant_qc.n_non_ref > 0)

    mtfails.write(mtout, overwrite = True)


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

    alspac_mt_giab_sample = annotate_mtfile_rf(mtfile, rf_htfile)
    alspac_mt_giab_sample_hard_filters = apply_hard_filters(alspac_mt_giab_sample)

    giab_vcf = resourcedir + "HG001_GRCh38_benchmark.interval.illumina.vcf.gz"
    giab_ht = prepare_giab_ht(giab_vcf, mtdir)

    mtout = mtdir + 'giab_vars_failing_hard_filters.mt'

    get_vars_in_giab_fail_hard_filters(alspac_mt_giab_sample_hard_filters, giab_ht, mtout, mtdir)


if __name__ == '__main__':
        main()