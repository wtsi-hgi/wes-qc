# useful functions
import datetime
import os
import sys
import yaml
import hail as hl
import pandas as pd
from shutil import rmtree
from typing import Optional, Union, Set
from gnomad.resources.resource_utils import TableResource


def str_timedelta(delta: datetime.timedelta) -> str:
    seconds_in_hr = 60 * 60
    seconds_in_day = seconds_in_hr * 24
    seconds = int(delta.total_seconds())
    days = seconds // seconds_in_day
    hours = (seconds % seconds_in_day) / seconds_in_hr
    return f"{days} days and {hours:4.2f} hr"


def get_script_path():
    #returns the path of the script that is being run
    return os.path.dirname(os.path.realpath(sys.argv[0]))

def parse_config():
    script_dir = get_script_path()
    input_yaml = script_dir + '/../config/inputs.yaml'
    if not os.path.exists(input_yaml):
        input_yaml = script_dir + '/../../config/inputs.yaml'
    with open(input_yaml, 'r') as y:
        inputs = yaml.load(y, Loader=yaml.FullLoader)

    return inputs


def expand_pd_array_col(
        df: pd.DataFrame,
        array_col: str,
        num_out_cols: int = 0,
        out_cols_prefix=None,
        out_1based_indexing: bool = True
) -> pd.DataFrame:
    """
    Expands a Dataframe column containing an array into multiple columns.
    :param DataFrame df: input dataframe
    :param str array_col: Column containing the array
    :param int num_out_cols: Number of output columns. If set, only the `n_out_cols` first elements of the array column are output.
                             If <1, the number of output columns is equal to the length of the shortest array in `array_col`
    :param out_cols_prefix: Prefix for the output columns (uses `array_col` as the prefix unless set)
    :param bool out_1based_indexing: If set, the output column names indexes start at 1. Otherwise they start at 0.
    :return: dataframe with expanded columns
    :rtype: DataFrame
    """

    if out_cols_prefix is None:
        out_cols_prefix = array_col

    if num_out_cols < 1:
        num_out_cols = min([len(x) for x in df[array_col].values.tolist()])

    cols = ['{}{}'.format(out_cols_prefix, i + out_1based_indexing) for i in range(num_out_cols)]
    df[cols] = pd.DataFrame(df[array_col].values.tolist())[list(range(num_out_cols))]

    return df


def get_rf(
    rf_dir: str,
    data: str = "rf_result",
    run_hash: Optional[str] = None,
) -> Union[str, TableResource]:
    """
    Gets the path to the desired RF data.
    Data can take the following values:
        - 'training': path to the training data for a given run
        - 'model': path to pyspark pipeline RF model
        - 'rf_result' (default): path to HT containing result of RF filtering
    :param str data: One of 'training', 'model' or 'rf_result' (default)
    :param str run_hash: Hash of RF run to load
    :return: Path to desired RF data
    """
    hashdir = rf_dir + run_hash + "/"
    model_file = hashdir + "rf.model"
    data_file = hashdir + data + ".ht"
    if data == "model":
        return model_file
    else:
        return TableResource(data_file)


def rm_mt(path: str):
    rmtree(path.replace('file:/', '/'))


def collect_pedigree_samples(ped: hl.Pedigree) -> Set[str]:
    samples = {getattr(trio, member) for trio in ped.trios for member in ('mat_id', 'pat_id', 's')}
    samples.discard(None)
    return samples

def clear_temp_folder(tmp_dir: str) -> None:
    if not tmp_dir.startswith("file://"):
        return
    tmp_dir = tmp_dir.replace("file://", "")
    print(f"=== Cleaning up temporary folder {tmp_dir}")
    try:
        shutil.rmtree(tmp_dir)
    except FileNotFoundError:
        print("TMP folder was cleaned up by Spark. Proceeding to the next step.")
    os.makedirs(tmp_dir, exist_ok=True)


def trace_analysis_step(mt: hl.MatrixTable, description: str = "") -> hl.MatrixTable:
    """
    The function records the name of the caller function into the global annotation of the matrixtable
    """
    caller_name = currentFuncName(1)
    step_description = hl.struct(step_name=caller_name, step_description=description)
    print(f"=== Annotating step {caller_name}")
    if "analysis_steps" in mt.globals:
        steps = mt.globals.analysis_steps
        steps = steps.append(step_description)
    else:
        steps = hl.array([step_description])
    print_analysis_steps(steps)
    mt = mt.annotate_globals(analysis_steps=steps)
    return mt

def select_founders(ped: hl.Pedigree) -> Set[str]:
    samples = collect_pedigree_samples(ped)
    for trio in ped.trios:
        samples.discard(trio.s)
    return samples
