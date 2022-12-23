# Filter out samples which have failed sample QC
# Filter to remove samples sequenced at Broad
import hail as hl
import pyspark
from wes_qc.utils.utils import parse_config


def filter_to_sanger_only(annotated_mt_file: str, sanger_mt_file: str):
    '''
    param str annotated_mt_file: File containing MatrixTable annotated with sequencing location
    param str sanger_mt_file: File for sanger only MatrixTable
    population and run id but not filtered
    :return: MatrixTable containing only Sanger samples
    :rtype: hl.MatrixTable
    '''
    mt = hl.read_matrix_table(annotated_mt_file)
    mt = mt.filter_cols(mt.sequencing_location == 'Sanger')  # filter to Sanger only
    mt.write(sanger_mt_file)


def remove_sample_qc_fails(sanger_mt_file: str, qc_filter_ht_file: str, samples_failing_qc_file: str, sample_qc_filtered_mt_file: str):
    '''
    param str sanger_mt_file: Input MatrixTable file
    param str qc_filter_ht_file: File contaning sample QC output
    param str samples_failing_qc_file: Output file for list of samples failing QC
    param str sample_qc_filtered_mt_file: Output file for MatrixTable with sample QC fails removed
    '''
    sangermt = hl.read_matrix_table(sanger_mt_file)
    sample_qc_ht = hl.read_table(qc_filter_ht_file)
    # identify samples which have failed any of the metrics tested
    sample_qc_ht = sample_qc_ht.annotate(filter_fail_count=(hl.len(sample_qc_ht.qc_metrics_filters)))

    filter_expr = (hl.case()
                   .when(sample_qc_ht.filter_fail_count == 0, 'Pass')
                   .when(sample_qc_ht.filter_fail_count > 0, 'Fail')
                   .default("")
                   )
    sample_qc_ht = sample_qc_ht.annotate(filter_result=filter_expr).key_by('s')
    # save a list of samples that have failed QC
    fail_ht = sample_qc_ht.filter(sample_qc_ht.filter_result == 'Fail').key_by()
    fails = fail_ht.select(fail_ht.s)
    fails.export(samples_failing_qc_file)
    # filter the sangermt to remove samples that have failed
    fails = fails.key_by('s')
    sangermt = sangermt.filter_cols(hl.is_defined(fails[sangermt.s]), keep=False)
    # save the filtered mt to file
    sangermt.write(sample_qc_filtered_mt_file, overwrite=True)


def main():
    # set up
    inputs = parse_config()
    mtdir = inputs['matrixtables_lustre_dir']
    annotdir = inputs['annotation_lustre_dir']

    # initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    qc_filter_ht_file = mtdir + "mt_pops_QC_filters.ht"
    annotated_mt_file = mtdir + "dv_unprocessed_with_pop.mt"  # annotated but unfiltered mt
    sample_qc_filtered_mt_file = mtdir + "mt_pops_QC_filters_sequencing_location_and_superpop_sanger_only_after_sample_qc.mt"
    sanger_mt_file = mtdir + "mt_pops_QC_filters_sequencing_location_and_superpop_sanger_only_after_sample_qc_sanger_only.mt"
    filter_to_sanger_only(annotated_mt_file, sanger_mt_file)
    samples_failing_qc_file = annotdir + "samples_failing_qc.tsv.bgz"
    remove_sample_qc_fails(sanger_mt_file, qc_filter_ht_file, samples_failing_qc_file, sample_qc_filtered_mt_file)


if __name__ == '__main__':
    main()
