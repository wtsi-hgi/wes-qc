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

def get_script_path() -> str:
    #returns the path of the script that is being run
    return os.path.dirname(os.path.realpath(sys.argv[0]))

# TODO: move this stuff to another utils file. Something like `common.py`

def getp(_dict: dict, keypath: str, silent: bool = False, default = None):
    """
    getPath: Get a value from a tree of dicts using dot-notation. That is, 
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

def checkp(_dict: dict, keypath: str):
    """
    checkPath: Check if a value is in a nested dict using dot-notation. That is, 
    `isinp(conf, "a.b.c")` is True if `conf['a']['b']['c']` is valid.
    That is, `conf['a']` and `conf['a']['b']` are dicts and `conf['a']['b']['c']` is 
    either a dict or a value.  
    This function will not check if a value is a dict or a field.

    Parameters
    - `_dict` is a nested dict
    - `keypath` is a sequence of keys joined by '.'
    """
    keyseq = keypath.strip('.').split('.')
    if keyseq[0] not in _dict:
        return False
    
    section = _dict
    breadcrumbs = ''
    for nextkey in keyseq:
        if isinstance(section, dict):
            if nextkey in section:
                section = section[nextkey]
                breadcrumbs += '.' + nextkey
            else:
                return False
        else:
            return True
    return True

def subdict(_dict: dict, keys: Iterable[str]) -> dict:
    """
    Get a subset of a dictionary safely. Output dictionary will only have those keys that are in both _dict and keys list.

    Useful for constructing a kwargs dict
    """
    correct_keys = _dict.keys() & set(keys)
    return {k: _dict[k] for k in correct_keys}

def multigetp(_dict: dict, keypaths: Iterable[str], silent: bool = False, default = None) -> tuple:
    """
    MultiGetPath: Get a series of items by keys using `utils.getp`.
    Useful when you want to get several items at once (via unpacking).
    """
    return (getp(_dict, k, silent, default) for k in keypaths)

import collections
def flatten(dictionary: collections.abc.MutableMapping, parent_key=False, separator='.') -> dict:
    """
    Turn a nested dictionary into a flattened dictionary.

    Maintains an insertion order.
    From https://stackoverflow.com/a/62186053

    :param dictionary: The dictionary to flatten
    :param parent_key: The string to prepend to dictionary's keys.
    :param separator: The string used to separate flattened keys
    :return: A flattened dictionary
    """

    items = []
    for key, value in dictionary.items():
        new_key = str(parent_key) + separator + key if parent_key else key
        if isinstance(value, collections.abc.MutableMapping):
            items.extend(flatten(value, new_key, separator).items())
        else:
            items.append((new_key, value))
    return dict(items)

def flat_to_nested(flat_dictionary: collections.abc.MutableMapping, parent_key=None, separator='.') -> dict:
    """
    Turn a flat dictionary with `separator`-separated keys into a nested dictionary.
    Maintains an insertion order only for fields.

    :param dictionary: A flat dictionary to transform to a nested one 
    :param parent_key: The string to prepend to dictionary's keys
    :param separator: The string used to separate flattened keys (default ".", a dot)
    :return: A nested dictionary
    """

    nested_dict = dict()

    for key, value in flat_dictionary.items():
        keyp = key.strip(separator).split(separator)
        section = nested_dict
        for k in keyp[:-1]:
            if k not in section:
                section[k] = dict()
            section = section[k]
        section[keyp[-1]] = value
        
    return nested_dict

def is_subsection(subsection_key: str, supersection_key: str) -> bool:
    """
    Compare two string keys for a nested directory in a dot-notation.
    Return True if `subsection_key` is a possible field of a dict at `supersection_key`.

    Examples: 
    - `is_subsection('foo.bar.baz', 'foo.bar') == True`
    - `is_subsection('foo.bar.baz', 'foo') == True`
    - `is_subsection('any.other.key', '') == True`
    """
    if supersection_key == '': 
        return True
    return subsection_key.startswith(supersection_key + '.')

def parent_section(section_key: str) -> str:
    """
    Extract a parent section for a nested dictionary string key in a dot notation.

    Examples: 
    - `parent_section('foo.bar.baz') == 'foo.bar'`
    - `parent_section('foo') == ''`
    - `parent_section('') == ''`
    """
    dot_idx = section_key.rfind('.')
    if dot_idx < 0: return ''
    return section_key[:dot_idx]

def resolve_cvar(cvar: str, base_config: dict, fieldname: str = '', undefined_ok: bool = False) -> str:
    """
    Get a value that cvar (without braces) should be substituted with.
    
    Parameters:
    - cvar: a cvar expression without braces. A field key in a nested dict, or a cvar constant
    - base_config: a context for cvar resolution. cvar will be resolved using data from base_config
    - fieldname: an optional key of the field that contains this cvar (in a dot-notation).  
        This is used to resolve local config variables and for logging
    - undefined_ok: If True and a cvar expression cannot be resolved, then return it as it is.  
        If False, throw a ValueError in that case.

    - If cvar is a cvar shortcut (a name defined in cvars section), then use the field it refers to
    - If cvar is a field name then 
        - try to use a field with that key in the same section (you can skip this by starting your cvar with a dot "{.foo.bar}")
        - try to find that field in the whole config 
    - If cannot find a suitable value and undefined_ok is False, raise a ValueError
    - If cannot find a suitable value and undefined_ok is True, return None 
    """
    try:    
        # 1. cvar shortcuts. Always absolute. At most once.
        cvar_shortcut_def = f"cvar.{cvar}"
        
        if checkp(base_config, cvar_shortcut_def):
            if undefined_ok:
                return getp(base_config, getp(base_config, cvar_shortcut_def), True, '{'+cvar+'}')
            return getp(base_config, getp(base_config, cvar_shortcut_def))

        # 2. config variable in the same section, providing we have a fieldname.
        #    skip if cvar is explicit full key (starts with `.`) or field is top-level (no `.`).
        if cvar[0] != '.' and '.' in fieldname:
            cvar_as_local = f"{parent_section(fieldname)}.{cvar}"
            if fieldname and (cvar_as_local in base_config):
                return getp(base_config, cvar_as_local)
            # note that we still need to check a full key
        
        # 3. config variable with a full key
        if cvar[0] == '.':
            cvar = cvar[1:]
        if undefined_ok:
            return getp(base_config, cvar, True, '{'+cvar+'}')
        return getp(base_config, cvar)

    except KeyError as e:
        raise ValueError(f"cannot substitute undefined config variable '{cvar}' referenced in field '{fieldname}'")

"""
Detect field names that end in dir, file, indir, outdir, infile, outfile + optionally _local.
Capture groups are (in/out/None, dir/file, _local/None).  
The first group does not matter for a match, but is useful for parsing. 
"""
__is_path_field_re = re.compile(r"(out|in)?(dir|file)(_local)?$")
def __is_path_field(fieldname: str):
    return __is_path_field_re.search(fieldname) is not None

def _expand_cvars_str(base_config: dict, str_with_cvar: str, fieldname: str, as_path: bool = False, undefined_ok: bool = False):
    """
    Substitute config variables in a string.  
    Will follow cvar_shortcuts, this helps avoid specifying the full field name for frequently used fields.  
    If undefined_ok, ignore errors, leave invalid cvars as is.

    :param dict base_config: a flat config dictionary and a source for values to substitute
    :param str str_with_cvar: a string with a variable that will be substituted by value from `base_config`
    :param str fieldname: a name of this string field for logging purposes
    :param bool as_path: shall this string be treated as a file/dir path
    :param bool undefined_ok: shall we ignore undefined variables, or raise an exception
    
    Perform basic path normalization if `as_path`.
    
    Will use `cvars` section from the `base_config` first, 
    and literal field names from the `base_config` second.
    """

    # Find all cvars that look like {something}
    # - no dots at the begininng or at the end
    # - no more 1 dot in a row
    # - a-zA-Z0-9_ as identifiers
    cvar_re = re.compile(r"\{([a-zA-Z0-9_]+(?:\.[a-zA-Z0-9_]+)*)\}")
    def repl(cvar_match):
        cvar_expr = cvar_match.group(1)
        return str(resolve_cvar(cvar_expr, base_config, fieldname, undefined_ok))
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

#TODO mk43: MERGE
def _process_cvars_in_config(unprocessed_config: dict, base_config: dict = dict()):
    """
    Recursively expand cvars in all string fields in a FLAT `unprocessed_config`
    using data from `base_config`.  
    Detect and activate path mode automatically based on __is_path_field.
    Specify additional config field aliases in `custom_shortcuts`.

    :param dict unprocessed_config: a flat dictionary with string values that needs cvar expansion
    :param dict base_config: a read-only dictionary, a base config, that will act as an additional source for config variables
    :return dict: 

    Order of config variable resolution:
    1. `cvar_shortcuts` section from the unprocessed_config
    2. `cvar_shortcuts` section from the base_config
    3. unprocessed_config fields in same section
    4. base_config fields in same section
    5. unprocessed_config fields by full key
    6. base_config fields by full key
    """

    raise NotImplementedError()

    return config

#TODO mk43: MERGE
def parse_config_file(file_obj, additional_cvar_shortcuts: dict = dict()):
    config = yaml.safe_load(file_obj)
    if 'cvar_shortcuts' not in config:
        config['cvar_shortcuts'] = dict()
    config['cvar_shortcuts'].update(additional_cvar_shortcuts)
    config = flatten(config)
    config = _process_cvars_in_config(config)
    return config

#TODO mk43: MERGE
# TODO: rename to load_config or get_config
def parse_config(path: str = None, additional_cvar_shortcuts: dict = dict()):
    """
    Read and parse WES QC YAML config from default location, env variable, or function argument.  
    Does variable substitution via `_process_cvars_in_config`.
    
    Config lookup order:
    1. `path` if it is not None
    2. `$WES_CONFIG` if `WES_CONFIG` env variable starts with `/` or `./`
    3. `../config/$WES_CONFIG` and `../../config/$WES_CONFIG`
    4. `../config/input.yaml` and `../../config/input.yaml`
    """

    # with an option to get the config file path from env or from arg
    
    if path is not None:
        print(f"Loading config '{path}', function arg")
        with open(path, 'r') as y:
            config = parse_config_file(y, additional_cvar_shortcuts)
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
        config = parse_config_file(y, additional_cvar_shortcuts)

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
