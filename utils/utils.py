# useful functions
import os
import sys
import hail as hl
from gnomad.resources.resource_utils import TableResource
import pandas as pd
from shutil import rmtree
from typing import Optional, Union, Set

"""
Config
"""

# TODO: cleanup these imports
from utils.config import (
    default_cvars,
    getp,
    subdict,
    multigetp,
    __is_path_field_re,
    __is_path_field,
    _expand_cvars,
    _expand_cvars_recursively,
    parse_config,
    path_local,
    path_spark
)

"""
Other
"""

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

import subprocess

# === Utils for downloading test data from the s3 storage === #

TEST_DATA_DOWNLOAD_URL = 'https://wes-qc-data.cog.sanger.ac.uk/all_test_data/test_data.zip'

# TODO: make versatile, don't download if already exists
def download_test_data_from_s3(outdir: str, move_dirs: dict, clean_up_unzip_dir: bool = False) -> None:
    """Download compressed test data from the s3 storage and 
    move the subfolders to the correct destinations inside the repo.

    Parameters
    ----------
    outdir : str
        Path to download and extract compressed test data.
    
    move_dirs : dict
        Directories to move. Keys are the dirs to be moved, values are their destination dirs
    
    clean_up_unzip_dir : bool, default: False
        Remove the folder with the unzipped data after moving the data from it.

    Return
    ------
        None
    """
    print(f'Downloading data from the s3 bucket') # TODO: improve logging
    # Download zipped archive with all the data
    subprocess.run(['wget', '-nc', TEST_DATA_DOWNLOAD_URL, '-P', outdir]) # skip existing

    # Unzip the data
    print(f'Unzipping the data')
    subprocess.run(['unzip', '-n', os.path.join(outdir, 'test_data.zip'), '-d', outdir])
    # Move unzipped data to correct folders in the test dir
    print(f'Moving data to correct dirs')

    for dir_to_move, destination_dir in move_dirs.items():
        subprocess.run(['mv', '-vn', dir_to_move, destination_dir])

    if clean_up_unzip_dir:
        subprocess.run(['rm', '-r', outdir])
