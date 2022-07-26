#get medians for number of variants per sample at each rf bin and for each cq
import hail as hl
import pyspark
import argparse
from wes_qc.utils.utils import parse_config
import datetime


def get_options():
    '''
    Get options from the command line
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("--runhash", help="RF training run hash")
    args = parser.parse_args()
    if not args.runhash:
        print("--runhash must be specified")
        exit(1)

    return args


def get_median_vars_per_sample(mtfile: str, bin: int, consequence: str) -> float:
    '''
    get median numbers of variants per sample for any given consequence and random forest bin
    :param str mtfile: random forest and consequence annotated mtfile
    :param int bin: maximum bin number
    :param str consequence: vep consequence
    :return: float median variants per sample with the specified consequence and maximum bin
    '''
    mt = hl.read_matrix_table(mtfile)
    mt = mt.filter_rows((mt.info.consequence == consequence) & (mt.info.rf_bin <= bin))
    mt = hl.sample_qc(mt)
    sampleqc_ht = mt.cols()
    med_vars_per_sample = sampleqc_ht.aggregate(hl.agg.approx_quantiles(sampleqc_ht.sample_qc.n_non_ref, 0.5))

    return med_vars_per_sample


def main():
    # set up
    args = get_options()
    inputs = parse_config()
    mtdir = inputs['matrixtables_lustre_dir']

    # initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    mtfile = mtdir + "mt_varqc_splitmulti_with_cq_and_rf_scores_" + args.runhash + ".mt"
    bin = 40
    consequence = 'missense_variant'
    print(datetime.datetime.now())
    median_bin_cq = get_median_vars_per_sample(mtfile, bin, consequence)
    print(isinstance(median_bin_cq, float))
    print(median_bin_cq)
    print(datetime.datetime.now())

if __name__ == '__main__':
    main()
