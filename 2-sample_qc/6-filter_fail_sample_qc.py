# Filter out samples which have failed sample QC
# Filter to remove samples sequenced at Broad
import hail as hl
import pyspark
import argparse
from utils.utils import parse_config

def get_options():
    '''
    Get options from the command line
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--add",
        help="add control samples", action="store_true")
    parser.add_argument("-f", "--filter",
        help="filter out a list of assigned poppulations", action="store_true")

    args = parser.parse_args()

    return args

def get_list_from_file(input_f: str):
        with open(input_f, 'r') as file:
            pop_list = [line.strip() for line in file]
        return pop_list

def filter_to_sas_only(mt: hl.MatrixTable, sp_list):
    for element in sp_list:
        mt = mt.filter_cols(mt.assigned_pop != element)
    #mt = mt.filter_cols(mt.assigned_pop != 'other-sas')#It was decided to remove only non-sas
    return mt


def remove_sample_qc_fails(sanger_mt_file: str, qc_filter_ht_file: str, samples_failing_qc_file: str, sample_qc_filtered_mt_file: str, cs_mt_file: str, add_control: bool, sp_to_remove):
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
    if len(sp_to_remove)>0:
        sangermt = filter_to_sas_only(sangermt, sp_to_remove)
    if add_control:
        cs_mt = hl.read_matrix_table(cs_mt_file)
        cs_mt = cs_mt.annotate_cols(assigned_pop='control')
        sangermt=sangermt.union_cols(cs_mt)
    # save the filtered mt to file
    sangermt.write(sample_qc_filtered_mt_file, overwrite=True)


def main():
    # set up
    args = get_options()
    inputs = parse_config()
    mtdir = inputs['matrixtables_lustre_dir']
    annotdir = inputs['annotation_lustre_dir']

    if args.filter:
        sp_file=inputs['populations_to_remove']
        sp_to_remove=get_list_from_file(sp_file)
    else:
        sp_to_remove=[]
    # initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")
    cs_mt_file = mtdir + "gatk_unprocessed.control_samples.mt"
    qc_filter_ht_file = mtdir + "mt_pops_QC_filters.ht"
    annotated_mt_file = mtdir + "gatk_unprocessed_with_pop.mt"  # annotated but unfiltered mt
    sample_qc_filtered_mt_file = mtdir + "mt_pops_QC_filters_sample_qc.mt"
    samples_failing_qc_file = annotdir + "samples_failing_qc.tsv"
    remove_sample_qc_fails(annotated_mt_file, qc_filter_ht_file, samples_failing_qc_file, sample_qc_filtered_mt_file, cs_mt_file, args.add, sp_to_remove)


if __name__ == '__main__':
    main()
