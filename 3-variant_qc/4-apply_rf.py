#apply RF model for variant QC
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
    parser.add_argument("--runhash", help="RF training run hash")
    args = parser.parse_args()
    if not args.runhash:
        print("--runhash must be specified")
        exit(1)

    return args


def main():
    # set up
    args = get_options()
    inputs = parse_config()
    rf_dir = inputs['var_qc_rf_dir']
    mtdir = inputs['matrixtables_lustre_dir']
    test_interval = inputs['rf_test_interval']

    # initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    run_hash = args.runhash
    rf_model = load_model(get_rf(rf_dir, data="model", run_hash=run_hash))
    ht = get_rf(rf_dir, data="training", run_hash=run_hash).ht()
    features = hl.eval(ht.features)
    ht = apply_rf_model(ht, rf_model, features, label=constants.LABEL_COL)
    ht = ht.annotate_globals(rf_hash=run_hash)
    ht = ht.checkpoint(get_rf(rf_dir, "rf_result_", run_hash=run_hash).path, overwrite=True)
    ht_summary = ht.group_by("tp", "fp", constants.TRAIN_COL, constants.LABEL_COL, constants.PREDICTION_COL).aggregate(n=hl.agg.count())
    ht_summary.show(n=20)
    

if __name__ == '__main__':
    main()