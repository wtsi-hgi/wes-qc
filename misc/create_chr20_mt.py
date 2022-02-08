# Create a subset of the unprocessed matrix table from 
# 1-import_data/2-import_gatk_vcfs_to_hail.py for use in testing
import hail as hl
import pyspark
import yaml
import os
import sys

def get_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))

def main():
    #set up paths and initialise hail
    script_dir = get_script_path()
    input_yaml = script_dir + '/../config/inputs.yaml'
    with open(input_yaml, 'r') as y:
        inputs = yaml.load(y, Loader=yaml.FullLoader)
    mtdir = inputs['matrixtables_lustre_dir']
    lustre_dir = inputs['hail_lustre_dir']
    sc = pyspark.SparkContext()
    hl.init(sc=sc, tmp_dir=lustre_dir, local_tmpdir=lustre_dir, default_reference="GRCh38")

    mt_in_file = mtdir + "/gatk_unprocessed.mt"
    mt_out_file = mtdir + "/tests/gatk_unprocessed_chr20.mt"
    print("Reading input matrix")
    mt = hl.read_matrix_table(mt_in_file)
    intervals = ['chr20']
    print("Subsetting to chr20")
    mt = hl.filter_intervals(mt, [hl.parse_locus_interval(x, reference_genome='GRCh38') for x in intervals])
    mt.write(mt_out_file, overwrite=True)

if __name__ == '__main__':
    main() 