#make RF model for variant QC
import hail as hl
import pyspark
import wes_qc.constants as constants
from wes_qc.utils.utils import parse_config, get_rf


def main():
    #set up
    inputs = parse_config()
    mtdir = inputs['matrixtables_lustre_dir']
    resourcedir = inputs['resource_dir']
    print(constants.TRUTH_DATA)


if __name__ == '__main__':
    main()