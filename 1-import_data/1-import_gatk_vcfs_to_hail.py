# Load GATK VCFs into hail and save as matrixtable
import hail as hl
import pyspark
import yaml
import os
import sys
import re
from wes_qc.utils.utils import parse_config

def load_vcfs_to_mt(indir, outdir, tmp_dir, header):
    '''
    load VCFs and save as hail mt
    '''
    objects = hl.utils.hadoop_ls(indir)
    
    # for some reason, paths prefix is `file:`, not a `file://`
    vcf_pattern = re.compile("file:.*vcf.b?gz")
    
    vcfs = [vcf["path"] for vcf in objects if vcf_pattern.match(vcf["path"])]
    #create and save MT

    print(f"info: Found {len(vcfs)} VCFs in {indir}")
    if header:
        print("info: Loading VCFs with header")
        mt = hl.import_vcf(vcfs, array_elements_required=False, force_bgz=True, header_file = header)
    else:
        print("info: Loading VCFs WITHOUT header")
        mt = hl.import_vcf(vcfs, array_elements_required=False, force_bgz=True)
        
    print("Saving as hail mt")
    mt_out_file = os.path.join(outdir, "gatk_unprocessed.mt")
    mt.write(mt_out_file, overwrite=True)


def main():
    #set up input variables
    inputs = parse_config()
    # dict.get returns None on KeyError
    vcf_header = inputs['step1_import'].get('gatk_vcf_header')
    import_vcf_dir = inputs['step1_import']['gatk_vcf_dir']
    mtdir = inputs['general']['matrixtables_lustre_dir']

    #initialise hail
    tmp_dir = inputs['general']['tmp_dir']
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
