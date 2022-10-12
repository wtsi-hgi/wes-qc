#get variants in GIAB missing from ALSPAC
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
    mt = mt.filter_entries(mt.hard_filters == 'Pass')

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

    tmpht = mtdir + "tmpg.ht"
    giab_vars = giab_vars.checkpoint(tmpht, overwrite = True)

    return giab_vars


def get_missing_vars(giab_ht: hl.Table, alspac_mt: hl.MatrixTable, outfile: str, mtdir: str):
    '''
    Get variants from GIAB which are not in ALSPAC
    :param hl.Table giab_ht: GIAB ht
    :param hl.MatrixTable alspac_mt: ALPSPAC matrixtable
    :param str outfile: output file path
    :param str mtdir: matrixtable directory
    '''
    #checkpoint files for speed
    tmpht = mtdir + "tmp3.ht"
    giab_ht = giab_ht.filter(hl.is_snp(giab_ht.alleles[0], giab_ht.alleles[1]))
    giab_ht = giab_ht.checkpoint(tmpht, overwrite = True)

    alspac_ht = alspac_mt.rows()
    alspac_ht = alspac_ht.filter(hl.is_snp(alspac_ht.alleles[0], alspac_ht.alleles[1]))
    tmpht2 = mtdir + "tmp4.ht"
    alspac_ht = alspac_ht.checkpoint(tmpht2, overwrite = True)

    giab_only = giab_ht.anti_join(alspac_ht)
    tmpht3 = mtdir + "tmp5.ht"
    giab_only = giab_only.checkpoint(tmpht3, overwrite = True)
    varcount = giab_only.count()
    print(varcount)

    loci = giab_only.locus.collect()
    refs = giab_only.alleles[0].collect()
    alts = giab_only.alleles[1].collect()

    with open(outfile, 'w') as o:
        for i in range(0, len(loci)):
            outline = ("\t").join([loci[i], refs[i], alts[i]])
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

    rf_htfile = rf_dir + "6617f838" + "/_gnomad_score_binning_tmp.ht"
    mtfile = mtdir + "mt_varqc_splitmulti.mt"

    alspac_mt_giab_sample = annotate_mtfile_rf(mtfile, rf_htfile)
    alspac_mt_giab_sample_hard_filters = apply_hard_filters(alspac_mt_giab_sample)

    giab_vcf = resourcedir + "HG001_GRCh38_benchmark.interval.illumina.vcf.gz"
    giab_ht = prepare_giab_ht(giab_vcf, mtdir)

    outfile_no_hard_filters = plot_dir + "/giab_vars_missing_from_alspac_unfiltered.txt"
    outfile_hard_filters = plot_dir + "/giab_vars_missing_from_alspac_hard_filters.txt"
    get_missing_vars(giab_ht, alspac_mt_giab_sample, outfile_no_hard_filters, mtdir)
    get_missing_vars(giab_ht, alspac_mt_giab_sample_hard_filters, outfile_hard_filters, mtdir)



if __name__ == '__main__':
        main()