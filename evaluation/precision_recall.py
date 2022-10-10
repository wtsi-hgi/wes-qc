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
    mt.count()
    exit(0)


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


if __name__ == '__main__':
    main()