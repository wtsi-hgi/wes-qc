# apply RF model for variant QC
import os

import hail as hl
import pyspark
import argparse
from utils import utils
import utils.constants as constants
from gnomad.variant_qc.random_forest import apply_rf_model, load_model

from wes_qc import hail_utils


def main() -> None:
    # set up
    inputs = utils.parse_config()

    data_root: str = inputs["data_root"]
    rf_dir = os.path.join(data_root, inputs["var_qc_rf_dir"])
    run_hash = inputs["runhash"]

    # initialise hail
    sc = hail_utils.init_hl(inputs["tmp_dir"])

    rf_model = load_model(utils.get_rf_model("file://" + rf_dir, run_hash=run_hash))
    ht = utils.get_rf_data("file://" + rf_dir, data="training", run_hash=run_hash).ht()
    features = hl.eval(ht.features)
    ht = apply_rf_model(ht, rf_model, features, label=constants.LABEL_COL)
    ht = ht.annotate_globals(rf_hash=run_hash)
    ht = ht.checkpoint(utils.get_rf_data("file://" + rf_dir, "rf_result", run_hash=run_hash).path, overwrite=True)
    ht_summary = ht.group_by("tp", "fp", constants.TRAIN_COL, constants.LABEL_COL, constants.PREDICTION_COL).aggregate(
        n=hl.agg.count()
    )
    ht_summary.show(n=20)

    hail_utils.stop_hl(sc)


if __name__ == "__main__":
    main()
