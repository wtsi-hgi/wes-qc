#get medians for number of variants per sample at each rf bin and for each cq
import hail as hl
import pyspark
import argparse
import matplotlib.pyplot as plt
import os
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


def get_median_vars_per_sample_per_bin_cq(mtfile: str, bins: list, consequences: list, plot_dir: str):
    '''
    get median numbers of variants per sample for any given consequence and random forest bin
    :param str mtfile: random forest and consequence annotated mtfile
    :param list bins: list of maximum bin numbers
    :param list consequences: vep consequences
    :param str plot_dir: output directory for plots
    '''
    mt = hl.read_matrix_table(mtfile)
    for consequence in consequences:
        median_variants_per_sample = []
        mt_cq = mt.filter_rows(mt.info.consequence == consequence)
        for bin in bins:
            #mt_tmp = mt.filter_rows((mt.info.consequence == consequence) & (mt.info.rf_bin <= bin))
            mt_tmp = mt_cq.filter_rows(mt_cq.info.rf_bin <= bin)
            mt_tmp = hl.sample_qc(mt_tmp)
            sampleqc_ht = mt_tmp.cols()
            med_vars_per_sample = sampleqc_ht.aggregate(hl.agg.approx_quantiles(sampleqc_ht.sample_qc.n_non_ref, 0.5))
            median_variants_per_sample.append(med_vars_per_sample)
        
        print(median_variants_per_sample)
        print(bins)
        plot_median_vars_per_cq(median_variants_per_sample, bins, consequence, plot_dir)


def plot_median_vars_per_cq(median_variants_per_sample: list, bins: list, consequence: str, plot_dir: str):
    '''
    plot median number of variants per sample for a specific consequence and range of bins
    :param list median_variants_per_sample: median variants per sample for this consequence
    :param list bins: list of maximum bin numbers
    :param list consequence: vep consequence
    :param str plot_dir: output directory for plots
    '''
    plotfile = plot_dir + "/" + consequence + ".png"
    plt.plot(bins, median_variants_per_sample)
    plt.savefig(plotfile)
    


def main():
    # set up
    args = get_options()
    inputs = parse_config()
    mtdir = inputs['matrixtables_lustre_dir']
    root_plot_dir = inputs['plots_dir_local']

    # initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    mtfile = mtdir + "mt_varqc_splitmulti_with_cq_and_rf_scores_" + args.runhash + ".mt"
    bins = list(range(5,101,10))
    consequences = ['missense_variant']
    plot_dir = root_plot_dir + "/variants_per_cq/" + args.runhash
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    print(datetime.datetime.now())
    get_median_vars_per_sample_per_bin_cq(mtfile, bins, consequences, plot_dir)
    print(datetime.datetime.now())


if __name__ == '__main__':
    main()
