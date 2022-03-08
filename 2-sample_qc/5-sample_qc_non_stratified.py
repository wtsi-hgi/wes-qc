# perform hail sample QC on raw matrixtable and annotate with superpops and fail count from stratified qc
import hail as hl
import pyspark
from wes_qc.utils.utils import parse_config

def run_sample_qc(raw_mt_file, pops_qc_filter_file, mt_sample_qc_file, ht_sample_qc_file):
    '''
    Run sample QC on raw variants
    :param str raw_mt_file: raw mt file 
    :param str pops_qc_filter_file: output from stratified QC
    :param str mt_sample_qc_file: output MT
    :param str ht_sample_qc_file: output HT, cols only
    '''
    mt = hl.read_matrix_table(raw_mt_file)
    mt = hl.sample_qc(mt)
    stratified_ht = hl.read_table(pops_qc_filter_file)
    # the qc_metrics_filters column lists the metrics each sample has failed on, 
    # but all fail on n_deletion and n_insertion so we want the number of items 
    # in this set minus 2 for metrics failed
    stratified_ht = stratified_ht.annotate(num_fails=(hl.len(stratified_ht.qc_metrics_filters) - 2))
    #now transfer num fails and superpopulation annotation to the MT
    mt = mt.annotate_cols(assigned_pop=stratified_ht[mt.s].assigned_pop)
    mt = mt.annotate_cols(filter_fail_count=stratified_ht[mt.s].num_fails)
    print("Writing annotated MT to file")
    mt.write(mt_sample_qc_file, overwrite=True)
    # also write just the cols as an HT
    mt.cols().write(ht_sample_qc_file,  overwrite=True)


def main():
    #set up
    inputs = parse_config()
    mtdir = inputs['matrixtables_lustre_dir']

    #initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    #run sample QC
    raw_mt_file = mtdir + "gatk_unprocessed.mt"
    mt_sample_qc_file = mtdir + "non_stratified_sample_qc.mt"
    pops_qc_filter_file = mtdir + "mt_pops_QC_filters.ht"
    ht_sample_qc_file = mtdir + "non_stratified_sample_qc_cols.ht"
    run_sample_qc(raw_mt_file, pops_qc_filter_file, mt_sample_qc_file, ht_sample_qc_file)


if __name__ == '__main__':
    main()