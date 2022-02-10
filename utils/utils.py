# useful functions
import os
import sys
import yaml

def get_script_path():
    #returns the path of the script that is being run
    return os.path.dirname(os.path.realpath(sys.argv[0]))

def parse_config():
    script_dir = get_script_path()
    input_yaml = script_dir + '/../config/inputs.yaml'
    with open(input_yaml, 'r') as y:
        inputs = yaml.load(y, Loader=yaml.FullLoader)

    return inputs