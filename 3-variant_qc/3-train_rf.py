# make RF model for variant QC
import os

import hail as hl
import pyspark
import uuid
import json
from pathlib import Path
from typing import Dict, List, Tuple, Any
import utils.constants as constants
from utils import utils
from gnomad.utils.file_utils import file_exists
from gnomad.variant_qc.pipeline import train_rf_model
from gnomad.variant_qc.random_forest import pretty_print_runs, save_model

from wes_qc import hail_utils

spark_local_message = """!!! WARNING !!!
The gnomAD fucntion train_rf_model() is a bit buggy
and can work incorrectly in the parallel SPARK environment.

If the run of the function will fail with some weird message
(no space left on device, wrong imports, etc),
try running model training on the master node only:

PYTHONPATH=$(pwd):$PYTHONPATH PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit --master local[*]  3-variant_qc/3-train_rf.py
"""


def get_rf_runs(rf_json_fp: str) -> Dict[Any, Any]:
    """
    Loads RF run data from JSON file.
    :param rf_json_fp: File path to rf json file.
    :return: Dictionary containing the content of the JSON file, or an empty dictionary if the file wasn't found.
    """
    if file_exists(rf_json_fp):
        print(f"=== Loading RF models data from {rf_json_fp}")
        with hl.hadoop_open("file://" + rf_json_fp) as f:
            return dict(json.load(f))
    else:
        print(f"=== File {rf_json_fp} could not be found. Returning empty RF run hash dict.")
        return {}


def train_rf(ht: hl.Table, test_intervals_str: str) -> Tuple[hl.Table, pyspark.ml.PipelineModel]:
    """
    Train RF model
    :param hl.Table ht: Hail table containing input data
    :param str test_intervals: Test intervals
    :return: Hail table and RF model
    """
    features = constants.FEATURES
    print(f"=== Test_intervals: {test_intervals_str}")

    fp_expr = ht.fail_hard_filters
    tp_expr = ht.omni | ht.mills | ht.kgp_phase1_hc | ht.hapmap
    ht = ht.annotate(tp=tp_expr, fp=fp_expr)

    if isinstance(test_intervals_str, str):
        test_intervals = [test_intervals_str]
        test_intervals = [hl.parse_locus_interval(x, reference_genome="GRCh38") for x in test_intervals]
        print("=== Resulting intervals")
        print(hl.eval(test_intervals))
    else:
        test_intervals = test_intervals_str
    ht = ht.persist()
    print("=== Start training the model")
    rf_ht, rf_model = train_rf_model(
        ht,
        rf_features=features,
        tp_expr=ht.tp,
        fp_expr=ht.fp,
        fp_to_tp=1.0,
        num_trees=500,
        max_depth=5,
        test_expr=hl.literal(test_intervals).any(lambda interval: interval.contains(ht.locus)),
    )
    # fp to tp = Ratio of FPs to TPs for training the RF model
    # num trees is number of trees in the model
    # max depth = maximum tree depth in model

    ht = ht.join(rf_ht, how="left")

    return ht, rf_model


def get_run_data(
    transmitted_singletons: bool,
    adj: bool,
    vqsr_training: bool,
    test_intervals: List[str],
    features_importance: Dict[str, float],
    test_results: List[hl.tstruct],
) -> Dict[str, Any]:
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
        transmitted_singletons = False

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


def main() -> None:
    # set up
    print(spark_local_message)

    inputs = utils.parse_config()

    data_root: str = inputs["data_root"]
    dataset_name: str = inputs["dataset_name"]
    mtdir: str = os.path.join(data_root, inputs["matrixtables_lustre_dir"])
    rf_dir = os.path.join(data_root, inputs["var_qc_rf_dir"])
    test_interval = inputs["rf_test_interval"]

    # initialise hail
    sc = hail_utils.init_hl(inputs["tmp_dir"])

    # hash for new RF
    runs_json = os.path.join(rf_dir, "rf_runs.json")
    rf_runs = get_rf_runs(runs_json)
    run_hash = str(uuid.uuid4())[:8]
    while run_hash in rf_runs:
        run_hash = str(uuid.uuid4())[:8]

    # train RF
    input_ht_file = os.path.join(mtdir, "ht_for_RF_by_variant_type_all_cols.ht")
    input_ht = hl.read_table("file://" + input_ht_file)

    ht_result, rf_model = train_rf(input_ht, test_interval)

    print("=== Writing out ht_training data")
    ht_result = ht_result.checkpoint(
        utils.get_rf_data("file://" + rf_dir, data="training", run_hash=run_hash).path, overwrite=True
    )

    rf_runs[run_hash] = get_run_data(
        vqsr_training=False,
        transmitted_singletons=True,
        test_intervals=test_interval,
        adj=True,
        features_importance=hl.eval(ht_result.features_importance),
        test_results=hl.eval(ht_result.test_results),
    )
    rf_runs[run_hash]["features_importance"] = dict(
        rf_runs[run_hash]["features_importance"]
    )  # convert featues importance frozendict to dict for dumping to json

    with hl.hadoop_open("file://" + runs_json, "w") as f:
        json.dump(rf_runs, f)
    pretty_print_runs(rf_runs)
    save_model(rf_model, utils.get_rf_model("file://" + rf_dir, run_hash=run_hash), overwrite=True)
    hail_utils.stop_hl(sc)


if __name__ == "__main__":
    main()
