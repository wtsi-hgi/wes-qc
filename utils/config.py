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

import collections
def flatten(dictionary: collections.abc.MutableMapping, parent_key=False, separator='.') -> dict:
    """
    From https://stackoverflow.com/a/62186053

    Turn a nested dictionary into a flattened dictionary
    :param dictionary: The dictionary to flatten
    :param parent_key: The string to prepend to dictionary's keys
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

def is_subsection(subsection_key: str, supersection_key: str) -> bool:
    if supersection_key == '': 
        return True
    return subsection_key.startswith(supersection_key + '.')

def parent_section(section_key: str) -> str:
    dot_idx = section_key.rfind('.')
    if dot_idx < 0: return ''
    return section_key[:dot_idx]

def resolve_cvar(cvar: str, base_config: dict, fieldname: str = '', undefined_ok: bool = False) -> str:
    """
    Get a value that cvar should be substituted with
    """
    try:    
        # 1. cvar shortcuts. Always absolute. At most once.
        cvar_shortcut_def = f"cvar_shortcuts.{cvar}"
        if cvar_shortcut_def in base_config:
            if undefined_ok:
                return base_config.get(base_config[cvar_shortcut_def], '{'+cvar+'}')
            return base_config[base_config[cvar_shortcut_def]]

        # 2. config variable in the same section, providing we have a fieldname.
        #    skip if cvar is explicit full key (starts with `.`) or field is top-level (no `.`).
        if cvar[0] != '.' and '.' in fieldname:
            cvar_as_local = f"{parent_section(fieldname)}.{cvar}"
            if fieldname and (cvar_as_local in base_config):
                return base_config[cvar_as_local]
            # note that we still need to check a full key
        
        # 3. config variable with a full key
        if cvar[0] == '.':
            cvar = cvar[1:]
        if undefined_ok:
            return base_config.get(cvar, '{'+cvar+'}')
        return base_config[cvar]

    except KeyError as e:
        raise ValueError(f"cannot substitute undefined config variable '{cvar}' referenced in field '{fieldname}'")

"""
Detect field names that end in dir, file, indir, outdir, infile, outfile + optionally _local.
Capture groups are (in/out/None, dir/file, _local/None) 
"""
__is_path_field_re = re.compile(r"(dir|file)(_local)?$")
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
    
    Will config `cvar_shortcuts` from the `base_config` first, 
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

    # the full context for variable substitution
    context = base_config.copy()
    # we only need cvar shortcuts for now
    for k,v in unprocessed_config.items():
        if k.startswith('cvar_shortcuts.'):
            context[k] = v
    # note that context does not contain unprocessed_config at start, but will contain in the end
    config = unprocessed_config.copy()

    for key in config:
        val = config[key]
        if isinstance(val, str):
            processed_str = _expand_cvars_str(context, val, key, as_path=__is_path_field(key))
            context[key] = processed_str
            config[key] = processed_str
        elif isinstance(val, dict):
            raise ValueError(f"This function accepts only a flat config. Nested dict detected at {key=}")

    return config

def parse_config_file(file_obj, additional_cvar_shortcuts: dict = dict()):
    config = yaml.safe_load(file_obj)
    if 'cvar_shortcuts' not in config:
        config['cvar_shortcuts'] = dict()
    config['cvar_shortcuts'].update(additional_cvar_shortcuts)
    config = flatten(config)
    config = _process_cvars_in_config(config)
    return config

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
