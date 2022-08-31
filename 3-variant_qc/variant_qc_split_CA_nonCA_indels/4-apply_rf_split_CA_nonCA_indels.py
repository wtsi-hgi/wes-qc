#apply RF model for variant QC
#modified for separate RFs for C>A SNPs, non-C>A SNPs and indels
import hail as hl
import pyspark
import argparse
from wes_qc.utils.utils import parse_config, get_rf
import wes_qc.utils.constants as constants
from gnomad.variant_qc.random_forest import apply_rf_model, load_model


def get_options():
    '''
    Get options from the command line
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("--runhash_CA", help="RF training run hash for C>A SNPs")
    parser.add_argument("--runhash_nonCA", help="RF training run hash for non-C>A SNPs")
    parser.add_argument("--runhash_indels", help="RF training run hash for indels")
    args = parser.parse_args()
    if not args.runhash_CA and args.runhash_nonCA and args.runhash_indels:
        print("--runhash_CA, runhash_nonCA and runhash_indels must be specified")
        exit(1)

    return args


def main():
    # set up
    args = get_options()
    inputs = parse_config()
    rf_dir = inputs['var_qc_rf_dir']

    # initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    run_hashes = [args.runhash_CA, args.runhash_nonCA, args.runhash_indels]
    for run_hash in run_hashes:
        rf_model = load_model(get_rf(rf_dir, data="model", run_hash=run_hash))
        ht = get_rf(rf_dir, data="training", run_hash=run_hash).ht()
        features = hl.eval(ht.features)
        ht = apply_rf_model(ht, rf_model, features, label=constants.LABEL_COL)
        ht = ht.annotate_globals(rf_hash=run_hash)
        ht = ht.checkpoint(get_rf(rf_dir, "rf_result", run_hash=run_hash).path, overwrite=True)
        ht_summary = ht.group_by("tp", "fp", constants.TRAIN_COL, constants.LABEL_COL, constants.PREDICTION_COL).aggregate(n=hl.agg.count())
        ht_summary.show(n=20)
    

if __name__ == '__main__':
    main()