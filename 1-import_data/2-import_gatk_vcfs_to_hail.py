from hail import Table
import hail as hl
import pyspark
import yaml
import os
import sys

def get_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))


def load_vcfs_to_mt(indir, outdir, header):
    '''
    load VCFs and save as hail mt
    '''
    objects = hl.utils.hadoop_ls(indir)
    vcfs = [vcf["path"] for vcf in objects if (vcf["path"].startswith("file") and vcf["path"].endswith("vcf.gz"))]
    print("Loading VCFs")
    mt = hl.import_vcf(vcfs, array_elements_required=False, force_bgz=True, header_file = header)
    print("Saving as hail mt")
    mt_out_file = outdir + "/gatk_unprocessed.mt"
    mt.write(mt_out_file, overwrite=True)


def main():
    #set up directories
    tmp_dir = "hdfs://spark-master:9820/"
    temp_dir = "file:///home/ubuntu/data/tmp"
    script_dir = get_script_path()
    input_yaml = script_dir + '/../config/inputs.yaml'
    with open(input_yaml, 'r') as y:
        inputs = yaml.load(y, Loader=yaml.FullLoader)

    vcf_header = inputs['gatk_vcf_header']
    import_vcf_dir = inputs['gatk_import_lustre_dir']
    mtdir = inputs['matrixtables_lustre_dir']
    lustre_dir = inputs['hail_lustre_dir']

    #initialise hail
    sc = pyspark.SparkContext()
    hl.init(sc=sc, tmp_dir=lustre_dir, local_tmpdir=lustre_dir, default_reference="GRCh38")

    #load VCFs
    load_vcfs_to_mt(import_vcf_dir, mtdir, vcf_header)


if __name__ == '__main__':
    main() 