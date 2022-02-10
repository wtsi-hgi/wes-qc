#apply gnomad's hard filters and impute sex
#input gatk_unprocessed.mt from step 1.1
import hail as hl
import pyspark
import yaml
from utils.utils import get_script_path

def main():
    script_dir = get_script_path()
    input_yaml = script_dir + '/../config/inputs.yaml'
    with open(input_yaml, 'r') as y:
        inputs = yaml.load(y, Loader=yaml.FullLoader)
    mtdir = inputs['matrixtables_lustre_dir']
    print(mtdir)

if __name__ == '__main__':
    main() 

