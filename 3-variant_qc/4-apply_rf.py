# apply RF model for variant QC
import hail as hl
from utils.utils import parse_config, get_rf, path_spark
import utils.constants as constants
from gnomad.variant_qc.random_forest import apply_rf_model, load_model
from wes_qc import hail_utils


def main():
    # = STEP SETUP = #
    config = parse_config()
    tmp_dir = config["general"]["tmp_dir"]

    # = STEP PARAMETERS = #
    model_id = config["general"]["rf_model_id"]

    # = STEP DEPENDENCIES = #
    rf_dir = path_spark(config["general"]["var_qc_rf_dir"])

    # = STEP OUTPUTS = #
    ## No outputs for this step

    # = STEP LOGIC = #
    _ = hail_utils.init_hl(tmp_dir)
    rf_model = load_model(get_rf(path_spark(rf_dir), data="model", model_id=model_id))
    ht = get_rf(path_spark(rf_dir), data="training", model_id=model_id).ht()
    features = hl.eval(ht.features)
    ht = apply_rf_model(ht, rf_model, features, label=constants.LABEL_COL)
    ht = ht.annotate_globals(rf_hash=model_id)
    ht = ht.checkpoint(get_rf(path_spark(rf_dir), "rf_result", model_id=model_id).path, overwrite=True)
    ht_summary = ht.group_by("tp", "fp", constants.TRAIN_COL, constants.LABEL_COL, constants.PREDICTION_COL).aggregate(
        n=hl.agg.count()
    )
    ht_summary.show(n=20)


if __name__ == "__main__":
    main()
