#make RF model for variant QC
import hail as hl
import pyspark
import uuid
import json
from typing import Dict, List, Tuple
import wes_qc.utils.constants as constants
from wes_qc.utils.utils import parse_config, get_rf
from gnomad.utils.file_utils import file_exists
from gnomad.variant_qc.pipeline import train_rf_model
from gnomad.variant_qc.random_forest import pretty_print_runs, save_model


def get_rf_runs(rf_json_fp: str) -> Dict:
    """
    Loads RF run data from JSON file.
    :param rf_json_fp: File path to rf json file.
    :return: Dictionary containing the content of the JSON file, or an empty dictionary if the file wasn't found.
    """
    if file_exists(rf_json_fp):
        with hl.hadoop_open(rf_json_fp) as f:
            return json.load(f)
    else:
        print(
            f"File {rf_json_fp} could not be found. Returning empty RF run hash dict."
        )
        return {}


def train_rf(ht: hl.Table, test_intervals: str) -> Tuple[hl.Table, pyspark.ml.PipelineModel]:
    '''
    Train RF model
    :param hl.Table ht: Hail table containing input data
    :param str test_intervals: Test intervals
    :return: Hail table and RF model
    '''
    features = constants.FEATURES
    print("test_intervals")
    print(test_intervals)

    fp_expr = ht.fail_hard_filters
    tp_expr = ht.omni | ht.mills | ht.kgp_phase1_hc | ht.hapmap
    ht = ht.annotate(tp=tp_expr, fp=fp_expr)

    if isinstance(test_intervals, str):
        test_intervals = [test_intervals]
        test_intervals = [
            hl.parse_locus_interval(x, reference_genome="GRCh38")
            for x in test_intervals
        ]
        print("Resulting intervals")
        print(hl.eval(test_intervals))

    ht=ht.persist()
   
    rf_ht, rf_model = train_rf_model(
        ht,
        rf_features=features,
        tp_expr=ht.tp,
        fp_expr=ht.fp,
        fp_to_tp=1.0,
        num_trees=500,
        max_depth=5,
        test_expr=hl.literal(test_intervals).any(
            lambda interval: interval.contains(ht.locus)),
    )
    #fp to tp = Ratio of FPs to TPs for training the RF model
    #num trees is number of trees in the model
    #max depth = maximum tree depth in model

    ht = ht.join(rf_ht, how="left")

    return ht, rf_model


def get_run_data(
    transmitted_singletons: bool,
    adj: bool,
    vqsr_training: bool,
    test_intervals: List[str],
    features_importance: Dict[str, float],
    test_results: List[hl.tstruct],


) -> Dict:
    """
    Creates a Dict containing information about the RF input arguments and feature importance
    :param bool transmitted_singletons: True if transmitted singletons were used in training
    :param bool adj: True if training variants were filtered by adj
    :param bool vqsr_training: True if VQSR training examples were used for RF training
    :param List of str test_intervals: Intervals withheld from training to be used in testing
    :param Dict of float keyed by str features_importance: Feature importance returned by the RF
    :param List of struct test_results: Accuracy results from applying RF model to the test intervals
    :return: Dict of RF information
    """
    if vqsr_training:
        transmitted_singletons = None

    run_data = {
        "input_args": {
            "transmitted_singletons": transmitted_singletons,
            "adj": adj,
            "vqsr_training": vqsr_training,
        },
        "features_importance": features_importance,
        "test_intervals": test_intervals,
    }

    if test_results is not None:
        tps = 0
        total = 0
        for row in test_results:
            values = list(row.values())
            # Note: values[0] is the TP/FP label and values[1] is the prediction
            if values[0] == values[1]:
                tps += values[2]
            total += values[2]
        run_data["test_results"] = [dict(x) for x in test_results]
        run_data["test_accuracy"] = tps / total

    return run_data


def main():
    # set up
    inputs = parse_config()
    print(inputs)
    exit(0)
    rf_dir = inputs['var_qc_rf_dir']
    mtdir = inputs['matrixtables_lustre_dir']
    test_interval = inputs['rf_test_interval']

    # initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    # hash for new RF
    run_hash = str(uuid.uuid4())[:8]
    rf_runs = get_rf_runs(rf_dir)
    while run_hash in rf_runs:
        run_hash = str(uuid.uuid4())[:8]
    

    # train RF
    input_ht = mtdir + "ht_for_RF_by_variant_type_all_cols.ht"
    runs_json = rf_dir + "rf_runs.json"
    ht_result, rf_model = train_rf(input_ht, test_interval)
    print("Writing out ht_training data")
    ht_result = ht_result.checkpoint(get_rf(data="training", run_hash=run_hash).path, overwrite=True)

    rf_runs[run_hash] = get_run_data(
        vqsr_training=False,
        transmitted_singletons=True,
        test_intervals=test_interval,
        adj=True,
        features_importance=hl.eval(ht_result.features_importance),
        test_results=hl.eval(ht_result.test_results),
    )

    with hl.hadoop_open(runs_json, "w") as f:
        json.dump(rf_runs, f)
        pretty_print_runs(rf_runs)
        save_model(
            rf_model, get_rf(data="model", run_hash=run_hash), overwrite=True)


if __name__ == '__main__':
    main()