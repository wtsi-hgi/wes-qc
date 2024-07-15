# Load GATK VCFs into hail and save as matrixtable
import pyspark
import yaml
import os
import sys
import re
import hail as hl
from utils.utils import parse_config

# DEBUG: for some reason, paths prefix is `file:`, not a `file://`
VCF_PATTERN = re.compile("file:.*vcf.b?gz")

def load_vcfs_to_mt(indir, outdir, tmp_dir, header):
    '''
    load VCFs and save as hail mt
    '''
    objects = hl.utils.hadoop_ls(indir)

    # get paths to vcf files
    vcfs = [vcf["path"] for vcf in objects if VCF_PATTERN.match(vcf["path"])]

    # create and save MT

    print(f"info: Found {len(vcfs)} VCFs in {indir}")
    if header:
        print("info: Loading VCFs with header")
        mt = hl.import_vcf(vcfs, array_elements_required=False, force_bgz=True, header_file=header)
    else:
        print("info: Loading VCFs WITHOUT header")
        mt = hl.import_vcf(vcfs, array_elements_required=False, force_bgz=True)
        
    print("Saving as hail mt")

    mt_out_file = os.path.join(outdir, "gatk_unprocessed.mt")
    mt.write(mt_out_file, overwrite=True)

def main():
    # set up input variables
    inputs = parse_config()
    # dict.get returns None on KeyError
    vcf_header = inputs.get('gatk_vcf_header')
    import_vcf_dir = inputs['gatk_import_lustre_dir']
    mtdir = inputs['matrixtables_lustre_dir']

    #initialise hail
    tmp_dir = inputs['tmp_dir']

    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    #load VCFs
    load_vcfs_to_mt(indir=import_vcf_dir, 
                    outdir=mtdir, 
                    tmp_dir=tmp_dir, 
                    header=vcf_header)

if __name__ == '__main__':
    main() 
