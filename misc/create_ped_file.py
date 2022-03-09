# Create ped file for trios to use in variant QC
from wes_qc.utils.utils import parse_config


def parse_manifest(manifest_file):
    '''
    Parse manifest file and return a dict keyed by proband.
    We only want complete trios.
    '''
    trios = {}
    with open(manifest_file, 'r') as f:
        lines = f.readlines()
        for l in lines:
            linedata = l.split("\t")
            print(linedata)
            exit(0)

    return trios


def write_ped():
    pass


def main():
    inputs = parse_config()
    resourcedir = inputs['resource_dir_local']
    manifest_file = resourcedir + "all_samples_with_proceed_and_seq_info_and_warehouse_info.txt"

    triodata = parse_manifest(manifest_file)
    write_ped()

if __name__ == '__main__':
    main()