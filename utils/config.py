"""
A separate file with minimal dependencies
"""

import os
import sys
import yaml
import re
from typing import Iterable

"""
Config utils
"""

def get_script_path():
    #returns the path of the script that is being run
    return os.path.dirname(os.path.realpath(sys.argv[0]))

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

# TODO: move this stuff to another utils file. Something like `common.py`
def subdict(_dict: dict, keys: Iterable[str]) -> dict:
    """
    Get a subset of a dictionary safely. Output dictionary will only have those keys that are in both _dict and keys list.

    Useful for constructing a kwargs dict
    """
    correct_keys = _dict.keys() & set(keys)
    return {k: _dict[k] for k in correct_keys}

def multigetp(_dict: dict, keypaths: Iterable[str], silent: bool = False, default = None) -> tuple:
    """
    Get a series of items by keys using `utils.getp`.
    Useful when you want to get several items at once (via unpacking).
    """
    return (getp(_dict, k, silent, default) for k in keypaths)

"""
Detect field names that end in dir, file, indir, outdir, infile, outfile + optionally _local.
Capture groups are (in/out/None, dir/file, _local/None) 
"""
__is_path_field_re = re.compile(r"(dir|file)(_local)?$")
def __is_path_field(fieldname: str):
    return __is_path_field_re.search(fieldname) is not None

def _expand_cvars(config: dict, str_with_cvar: str, fieldname: str, as_path: bool = False, custom_cvars: dict = None, undefined_ok: bool = False):
    """
    Expand config variables in a string.  
    Ignore errors, leave invalid cvars as is. 

    Perform basic path normalization if `as_path`.
    
    Will use `custom_cvars` first, config `cvars` second, and literal field names third.
    """
    # collect cvars from config
    cvars = config.get('cvars', dict())

    # Find all cvars that look like {something}
    # - no dots at the begininng or at the end
    # - no more 1 dot in a row
    # - a-zA-Z0-9_ as identifiers
    cvar_re = re.compile(r"\{([a-zA-Z0-9_]+(?:\.[a-zA-Z0-9_]+)*)\}")
    def repl(cvar_match):
        print(f"field {fieldname} {cvar_match.group(0)=}, ", end='')
        cvar_expr = cvar_match.group(1)
        if custom_cvars and cvar_expr in custom_cvars:
            print("custom cvar, ", end='')
            cvar_expr = custom_cvars[cvar_expr]
        elif cvar_expr in cvars: 
            print("cvar, ", end='')
            cvar_expr = cvars[cvar_expr]
        try:
            substitution = str(getp(config, keypath=cvar_expr, silent=undefined_ok, default='{'+cvar_expr+'}'))
            print(f'{cvar_expr} -> {substitution} ({undefined_ok=})')
        except KeyError as e:
            print(config)
            raise ValueError(f"field {fieldname}, expression {cvar_match.group(1)}, cannot substitute because variable is undefined.")
        return substitution
    expanded_str = cvar_re.sub(repl, str_with_cvar)

    if as_path:
        # separate protocol beforehand because os.path.normalize does bad things to URIs
        proto_re = re.compile(r"[a-zA-Z][a-zA-Z0-9]*:\/\/")
        proto_match = proto_re.match(expanded_str.strip())
        if proto_match is not None: 
            expanded_str = expanded_str[proto_match.end():]
        # normalize path
        expanded_str = os.path.normpath(expanded_str) if expanded_str else ''
        # restore the protocol
        if proto_match is not None:
            expanded_str = proto_match[0] + expanded_str

    return expanded_str

# TODO: fix infinite recursion
def _expand_cvars_recursively(config: dict, dict_to_expand, custom_cvars: dict = None, fieldname: str = ''):
    """
    Recursively expand cvars in all string fields in a nested dict structure `dict_to_expand`
    using data from `config`.  
    Detect and activate path mode automatically based on __is_path_field

    Will use `custom_cvars` first, `default_cvars` second, and literal field names third.
    """

    # A former version of this function had an in_place argument to perform this operation 
    # without creating a deep copy of the whole nested dictionary. This lead to errors in
    # variable substitution in cases when variable refers to a field that is defined 
    # further down the config file and requires a variable substitution on its own.
    # By constructing the parsed_config from the start via depth-first traversal we can 
    # detect such cases, as they lead to the KeyError in getp.
     
    parsed_config = dict()
    for key in config:
        next_fieldname = f"{fieldname}.{key}" if fieldname else key
        val = config[key]
        if isinstance(val, str):
            # substitute variables in the string
            parsed_config[key] = _expand_cvars(parsed_config, val, fieldname=next_fieldname, as_path=__is_path_field(key), custom_cvars=custom_cvars)
        elif isinstance(val, dict):
            # go deeper, substituting all variables in a subdict
            parsed_config[key] = _expand_cvars_recursively(parsed_config, val, custom_cvars=custom_cvars, fieldname=next_fieldname)
        else:
            # not a str and not a subdict, just copy the value
            parsed_config[key] = val
    return parsed_config

def parse_config(path: str = None, custom_cvars: dict = None):
    # with an option to get the config file path from env or from arg
    
    if path is not None:
        print(f"Loading config '{path}', function arg")
        with open(path, 'r') as y:
            inputs = yaml.load(y, Loader=yaml.FullLoader)
        config = _expand_cvars_recursively(inputs, inputs, custom_cvars=custom_cvars)
        return config
    
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
    config = _expand_cvars_recursively(inputs, inputs, inplace=True, custom_cvars=custom_cvars)

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
