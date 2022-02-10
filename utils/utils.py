# useful functions
import os
import sys

def get_script_path():
    #returns the path of the script that is being run
    return os.path.dirname(os.path.realpath(sys.argv[0]))
