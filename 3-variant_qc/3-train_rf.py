# make RF model for variant QC
import hail as hl
import pyspark
import uuid
import argparse
import json
from typing import Dict, List, Tuple
import utils.constants as constants
from utils.utils import parse_config, get_rf, path_spark
from gnomad.utils.file_utils import file_exists
from gnomad.variant_qc.pipeline import train_rf_model
from gnomad.variant_qc.random_forest import pretty_print_runs, save_model
from wes_qc import hail_utils

spark_local_message = """!!! WARNING !!!
The gnomAD fucntion train_rf_model() in some cases
could work incorrectly in the parallel SPARK environment.

If the run of the function will fail with some weird message
(no space left on device, wrong imports, etc),
try running model training on the master node only:

PYTHONPATH=$(pwd):$PYTHONPATH PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit --master local[*]  3-variant_qc/3-train_rf.py
"""


def get_options():
    """
    Get options from the command line
    """
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--manual-model-id",
        help="Manually set the identifier of the RF run. Use for reproducible runs. Don't specify to generate a new run hash.",
        type=str,
    )

    args = parser.parse_args()
    return args


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
        print(f"File {rf_json_fp} could not be found. Returning empty RF run hash dict.")
        return {}


def train_rf(ht: hl.Table, test_intervals: str, config: dict) -> Tuple[hl.Table, pyspark.ml.PipelineModel]:
    """
    Train RF model
    :param hl.Table ht: Hail table containing input data
    :param str test_intervals: Test intervals
    :return: Hail table and RF model
    """
    conf = config["step3"]["train_rf"]
    features = constants.FEATURES
    print("test_intervals")
    print(test_intervals)

    fp_expr = ht.fail_hard_filters
    tp_expr = ht.omni | ht.mills | ht.kgp_phase1_hc | ht.hapmap
    ht = ht.annotate(tp=tp_expr, fp=fp_expr)

    if isinstance(test_intervals, str):
        test_intervals = [test_intervals]
        test_intervals = [hl.parse_locus_interval(x, reference_genome="GRCh38") for x in test_intervals]
        print("Resulting intervals")
        print(hl.eval(test_intervals))

    ht = ht.persist()

    # TODO: use kwargs expantions from the config as this is a gnomad package function
    rf_ht, rf_model = train_rf_model(
        ht,
        rf_features=features,
        tp_expr=ht.tp,
        fp_expr=ht.fp,
        fp_to_tp=conf["gnomad_train_rf_fp_to_tp"],
        num_trees=conf["gnomad_train_rf_num_trees"],
        max_depth=conf["gnomad_train_rf_max_depth"],
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
    args = get_options()
    config = parse_config()
    tmp_dir = config["general"]["tmp_dir"]

    # = STEP PARAMETERS = #
    test_interval = config["step3"]["rf_test_interval"]  # used in multiple functions
    runs_json = config["step3"]["runs_json"]

    # = STEP DEPENDENCIES = #
    input_ht_file = config["step3"]["train_rf"]["input_ht_file"]

    # = STEP OUTPUTS = #
    rf_dir = path_spark(config["general"]["var_qc_rf_dir"])

    # = STEP LOGIC = #

    # initialise hail
    _ = hail_utils.init_hl(tmp_dir)

    # use manual hash value if provided to control the output folder name

    # new hash is generated by default
    model_id = str(uuid.uuid4())[:8]  # TODO: move to utils
    rf_runs = get_rf_runs(path_spark(rf_dir))
    if args.manual_model_id:
        print(f"Using provided run ID: {args.manual_model_id}")
        model_id = args.manual_model_id
    else:
        # hash for new RF
        while model_id in rf_runs:
            model_id = str(uuid.uuid4())[:8]

    # train RF
    input_ht = hl.read_table(path_spark(input_ht_file))

    ht_result, rf_model = train_rf(input_ht, test_interval, config)
    print("Writing out ht_training data")
    ht_result = ht_result.checkpoint(
        get_rf(path_spark(rf_dir), data="training", model_id=model_id).path, overwrite=True
    )

    rf_runs[model_id] = get_run_data(
        vqsr_training=False,
        transmitted_singletons=True,
        test_intervals=test_interval,
        adj=True,
        features_importance=hl.eval(ht_result.features_importance),
        test_results=hl.eval(ht_result.test_results),
    )
    rf_runs[model_id]["features_importance"] = dict(
        rf_runs[model_id]["features_importance"]
    )  # convert feature importance frozendict to dict for dumping to json

    print(f"Saving runs info in {runs_json}")
    # TODO: not sure if this works as intended. Shouldn't the run info json file be saved in the corresponding run folder?
    # If not, run info has to be appended, rather than overwriting existing runs.
    #
    with open(runs_json, "w") as f:
        json.dump(rf_runs, f)
    pretty_print_runs(rf_runs)

    save_model(rf_model, get_rf(path_spark(rf_dir), data="model", model_id=model_id), overwrite=True)


if __name__ == "__main__":
    print(spark_local_message)
    main()
