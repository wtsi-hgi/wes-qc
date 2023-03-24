# Filter out samples which have failed sample QC
# Filter to remove samples sequenced at Broad
import hail as hl
import pyspark
from utils.utils import parse_config


def filter_to_sas_only(mt: hl.MatrixTable):
    mt = mt.filter_cols(mt.assigned_pop != 'non-sas')
    mt = mt.filter_cols(mt.assigned_pop != 'other-sas')
    return mt


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
    fails = fail_ht.select(fail_ht.s, fail_ht.assigned_pop)
    print(f'{fails.count()} samples fail qc')
    fails.export(samples_failing_qc_file)
    # filter the sangermt to remove samples that have failed
    fails = fails.key_by('s')
    sangermt = sangermt.filter_cols(hl.is_defined(fails[sangermt.s]), keep=False)
    sangermt = filter_to_sas_only(sangermt)
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
    annotated_mt_file = mtdir + "gatk_unprocessed_with_pop.mt"  # annotated but unfiltered mt
    sample_qc_filtered_mt_file = mtdir + "mt_pops_QC_filters_sample_qc.mt"
    samples_failing_qc_file = annotdir + "samples_failing_qc.tsv"
    remove_sample_qc_fails(annotated_mt_file, qc_filter_ht_file, samples_failing_qc_file, sample_qc_filtered_mt_file)


if __name__ == '__main__':
    main()
