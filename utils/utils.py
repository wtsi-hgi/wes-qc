# useful functions
import os
import sys
import yaml
import hail as hl
import pandas as pd
from shutil import rmtree
from typing import Optional, Union, Set
from gnomad.resources.resource_utils import TableResource
import re

def get_script_path():
    #returns the path of the script that is being run
    return os.path.dirname(os.path.realpath(sys.argv[0]))


special_cvars = {
    'tmpdir': 'general.tmp_dir',
    'anndir': 'general.annotation_dir',
    'mtdir': 'general.matrixtables_dir',
    'resdir': 'general.resource_dir',
}


"""
Config utils
"""

def getp(_dict: dict, keypath: str, silent: bool = False, default = None):
    """
    Get a value from a tree of dicts using dot-notation. That is, 
    `getp(conf, "a.b.c") == conf['a']['b']['c']`

    - `_dict` is a structure of nested dicts
    - `keypath` is a sequence of keys joined by '.'
    - if `silent` then return `default` value on invalid path instead of raising an exception  

    Note that this function can return a dict as well as a scalar -- it has no output type checks.
    """
    keyseq = keypath.strip('.').split('.')
    if keyseq[0] not in _dict:
        if silent: return default
        raise KeyError(f"top-level dict has no '{keyseq[0]}' key")
    
    section = _dict
    breadcrumbs = ''
    for nextkey in keyseq:
        if isinstance(section, dict):
            if nextkey in section:
                section = section[nextkey]
                breadcrumbs += '.' + nextkey
            else:
                if silent: return default
                raise KeyError(f"{keypath} : '{breadcrumbs}' section has no '{nextkey}' field")
        else:
            if silent: return default
            raise KeyError(f"{keypath} : {breadcrumbs + '.' + nextkey} is a field, not a section")
    return section

"""
Detect field names that end in dir, file, indir, outdir, infile, outfile + optionally _local.
Capture groups are (in/out/None, dir/file, _local/None) 
"""
__is_path_field_re = re.compile(r"(out|in)?(dir|file)(_local)?$")
def __is_path_field(fieldname):
    return __is_path_field_re.search(fieldname) is not None

def __expand_cvars(config: dict, str_with_cvar: str, as_path: bool = False):
    """
    Expand config variables in a string.  
    Ignore errors, leave invalid cvars as is. 

    Perform basic path normalization if `as_path`
    """
    # Find all cvars that look like {cvar_name}
    # - no dots at the begininng or at the end
    # - no more 1 dot in a row
    # - a-zA-Z0-9_ as identifiers
    cvar_re = re.compile(r"\{([a-zA-Z0-9_]+(?:\.[a-zA-Z0-9_]+)*)\}")
    def repl(cvar_match):
        cvar_name = cvar_match.group(1)
        if cvar_name in special_cvars: 
            cvar_name = special_cvars[cvar_name]
        return str(getp(config, keypath=cvar_name, silent=True, default='{'+cvar_name+'}'))
    expanded_str = cvar_re.sub(repl, str_with_cvar)

    if as_path:
        # separate protocol beforehand because os.path.normalize does bad things to URIs
        proto_re = re.compile(r"[a-zA-Z][a-zA-Z0-9]*:\/\/")
        proto_match = proto_re.match(expanded_str.strip())
        if proto_match is not None: 
            expanded_str = expanded_str[proto_match.end():]
        # normalize path
        expanded_str = os.path.normpath(expanded_str)
        # restore the protocol
        if proto_match is not None:
            expanded_str = proto_match[0] + expanded_str

    return expanded_str

def __expand_cvars_recursively(config: dict, dict_to_expand, inplace=False):
    """
    Recursively expand cvars in all string fields in a nested dict structure `dict_to_expand`
    using data from `config`.  
    Detect and activate path mode automatically based on __is_path_field
    """
    if inplace:
        _dict = dict_to_expand
    else:
        from copy import deepcopy
        _dict = deepcopy(dict_to_expand)
    for key in _dict:
        val = _dict[key]
        if isinstance(val, str):
            _dict[key] = __expand_cvars(config, val, as_path=__is_path_field(key))
        elif isinstance(val, dict):
            # force inplace
            __expand_cvars_recursively(config, val, inplace=True)
    return _dict

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

    # expand config variables in strings
    config = __expand_cvars_recursively(inputs, inputs, inplace=True)

    return config


"""
Path utils
"""    

def path_spark(path: str):
    """
    Ensure `path` is a valid spark path. 
    That is, add 'file://' prefix if needed
    """
    if path.startswith('file://'):
        return path
    if path.startswith('hdfs://'):
        raise NotImplementedError('Cannot convert hdfs to a spark path')
    return 'file://' + path

def path_local(path: str):
    """
    Ensure `path` is a valid linux path. 
    That is, remove 'file://' prefix if needed
    """
    if path.startswith('file://'):
        return path[7:]
    if path.startswith('hdfs://'):
        raise NotImplementedError('Cannot convert hdfs to a local path')
    return path


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
