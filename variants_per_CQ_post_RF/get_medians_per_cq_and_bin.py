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


def get_median_vars_per_sample_per_bin_cq(mtfile: str, bins: list, consequences: list):
    '''
    get median numbers of variants per sample for any given consequence and random forest bin
    :param str mtfile: random forest and consequence annotated mtfile
    :param list bins: list of maximum bin numbers
    :param list consequences: vep consequences
    '''
    mt = hl.read_matrix_table(mtfile)
    median_variants_per_sample = []
    for consequence in consequences:
        for bin in bins:
            mt_tmp = mt.filter_rows((mt.info.consequence == consequence) & (mt.info.rf_bin <= bin))
            mt_tmp = hl.sample_qc(mt_tmp)
            sampleqc_ht = mt_tmp.cols()
            med_vars_per_sample = sampleqc_ht.aggregate(hl.agg.approx_quantiles(sampleqc_ht.sample_qc.n_non_ref, 0.5))
            median_variants_per_sample.append()

    print(bins)
    print(median_variants_per_sample)

    


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
    bins = list(range(1,101,5))
    consequences = ['missense_variant']
    print(datetime.datetime.now())
    get_median_vars_per_sample_per_bin_cq(mtfile, bins, consequences)
    print(datetime.datetime.now())


if __name__ == '__main__':
    main()
