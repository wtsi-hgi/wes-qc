# useful functions
import os
import sys
import yaml
import hail as hl
import pandas as pd
from shutil import rmtree
from typing import Optional, Union, Set
from gnomad.resources.resource_utils import TableResource

def get_script_path():
    #returns the path of the script that is being run
    return os.path.dirname(os.path.realpath(sys.argv[0]))

# TODO: caching
def parse_config():
    # with an option to get the config file path from env
    
    # find config dir
    script_dir = get_script_path()
    print("info: script_dir ", script_dir)
    config_dir = '../config'
    # to handle 3-variant_qc subdirs 
    if not os.path.exists(os.path.join(script_dir, config_dir)):
        config_dir = '../../config'
    config_dir = os.path.join(script_dir, config_dir)
    
    # default values
    input_yaml = os.path.join(config_dir, 'inputs.yaml')
    config_type = 'default'
    
    # handle environmental variable
    if 'WES_CONFIG' in os.environ:
        config_path = os.environ['WES_CONFIG']
        if config_path[0] == '/' or config_path.startswith('./'):
            # an absolute path
            input_yaml = config_path
            config_type = 'WES_CONFIG absolute'
        else:
            # consider path as relative to the config dir!
            input_yaml = os.path.join(config_dir, config_path)
            config_type = 'WES_CONFIG relative to config dir'

    # read config
    if not os.path.exists(input_yaml):
        print(f"error: config {input_yaml} does not exist")
    print(f"Loading config '{input_yaml}', {config_type}")
        
    with open(input_yaml, 'r') as y:
        inputs = yaml.load(y, Loader=yaml.FullLoader)

    return inputs

PATH_LOCAL, PATH_REMOTE, PATH_SAME = 1,2,4
def get_path(directory, filename, variant=PATH_SAME, user_prefix=None):
    dir_variant = PATH_REMOTE if directory.startswith('file://') else PATH_LOCAL
    if user_prefix:
        filename = f"{user_prefix}_" + filename
    if variant == PATH_SAME:
        pass
    if variant == PATH_REMOTE and dir_variant == PATH_LOCAL:
        directory = "file://" + directory
    if variant == PATH_LOCAL and dir_variant == PATH_REMOTE:
        directory = directory.replace("file://", "", count=1)
    return os.path.join(directory, filename)
        

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


def select_founders(ped: hl.Pedigree) -> Set[str]:
    samples = collect_pedigree_samples(ped)
    for trio in ped.trios:
        samples.discard(trio.s)
    return samples
