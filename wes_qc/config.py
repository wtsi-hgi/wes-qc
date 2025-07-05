import os
import sys
from typing import Optional, Any, cast
from omegaconf import OmegaConf

"""
Config utils - locate and render config files
"""

WxsqcConfig = dict[str, Any]


def get_script_path() -> str:
    # returns the path of the script that is being run
    return os.path.dirname(os.path.realpath(sys.argv[0]))


def parse_config_file(file_name: str) -> WxsqcConfig:
    config = OmegaConf.to_container(OmegaConf.load(file_name), resolve=True)
    return cast(WxsqcConfig, config)


def get_config(path: Optional[str] = None) -> WxsqcConfig:
    """
    Read and parse WES QC YAML config from the default location,
    env variable, or function argument.

    Config lookup order:
    1. `path` if it is not None
    2. `$WES_CONFIG` if `WES_CONFIG` env variable starts with `/` or `./`
    3. `../config/$WES_CONFIG` and `../../config/$WES_CONFIG`
    4. `../config/input.yaml` and `../../config/input.yaml`
    """

    # with an option to get the config file path from env or from arg

    if path is not None:
        print(f"Loading config '{path}', function arg")
        return parse_config_file(path)

    # find config dir
    script_dir = get_script_path()
    print("info: script_dir ", script_dir)
    config_dir = "../config"
    # to handle 3-variant_qc subdirs
    if not os.path.exists(os.path.join(script_dir, config_dir)):
        config_dir = "../../config"
    config_dir = os.path.join(script_dir, config_dir)

    # default values
    input_yaml = os.path.join(config_dir, "inputs.yaml")
    config_type = "default"

    # handle environmental variable
    if "WES_CONFIG" in os.environ:
        config_path = os.environ["WES_CONFIG"]
        if config_path[0] == "/" or config_path.startswith("./"):
            # an absolute path
            input_yaml = config_path
            config_type = "WES_CONFIG absolute"
        else:
            # consider path as relative to the config dir!
            input_yaml = os.path.join(config_dir, config_path)
            config_type = "WES_CONFIG relative to config dir"

    # read config
    if not os.path.exists(input_yaml):
        print(f"error: config {input_yaml} does not exist")
    print(f"Loading config '{input_yaml}', {config_type}")

    return parse_config_file(input_yaml)
