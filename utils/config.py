"""
A separate file with minimal dependencies
"""

import os
import yaml
import re

"""
A dictionary of predefined cvars in the format of  
`{ 'cvar_name': 'path.to.cvar.definition.in.config' }`
"""
default_cvars = {
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

def __expand_cvars(config: dict, str_with_cvar: str, as_path: bool = False, custom_cvars: dict = None):
    """
    Expand config variables in a string.  
    Ignore errors, leave invalid cvars as is. 

    Perform basic path normalization if `as_path`.
    
    Will use `custom_cvars` first, `default_cvars` second, and literal field names third.
    """
    # Find all cvars that look like {cvar_name}
    # - no dots at the begininng or at the end
    # - no more 1 dot in a row
    # - a-zA-Z0-9_ as identifiers
    cvar_re = re.compile(r"\{([a-zA-Z0-9_]+(?:\.[a-zA-Z0-9_]+)*)\}")
    def repl(cvar_match):
        cvar_name = cvar_match.group(1)
        if custom_cvars and cvar_name in custom_cvars:
            cvar_name = custom_cvars[cvar_name]
        elif cvar_name in default_cvars: 
            cvar_name = default_cvars[cvar_name]
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

def __expand_cvars_recursively(config: dict, dict_to_expand, inplace=False, custom_cvars: dict = None):
    """
    Recursively expand cvars in all string fields in a nested dict structure `dict_to_expand`
    using data from `config`.  
    Detect and activate path mode automatically based on __is_path_field

    Will use `custom_cvars` first, `default_cvars` second, and literal field names third.
    """
    if inplace:
        _dict = dict_to_expand
    else:
        from copy import deepcopy
        _dict = deepcopy(dict_to_expand)
    for key in _dict:
        val = _dict[key]
        if isinstance(val, str):
            _dict[key] = __expand_cvars(config, val, as_path=__is_path_field(key), custom_cvars=custom_cvars)
        elif isinstance(val, dict):
            # force inplace
            __expand_cvars_recursively(config, val, inplace=True, custom_cvars=custom_cvars)
    return _dict

def parse_config(path: str = None, custom_cvars: dict = None):
    # with an option to get the config file path from env or from arg
    
    if path is not None:
        print(f"Loading config '{path}', function arg")
        with open(path, 'r') as y:
            inputs = yaml.load(y, Loader=yaml.FullLoader)
        config = __expand_cvars_recursively(inputs, inputs, inplace=True, custom_cvars=custom_cvars)
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
    config = __expand_cvars_recursively(inputs, inputs, inplace=True, custom_cvars=custom_cvars)

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
