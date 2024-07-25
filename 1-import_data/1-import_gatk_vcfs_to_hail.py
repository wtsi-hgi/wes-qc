# Load GATK VCFs into hail and save as matrixtable
import hail as hl
import pyspark
import yaml
import os
import sys
import re
from wes_qc.utils.utils import parse_config, path_local, path_spark

def load_vcfs_to_mt(config):
    '''
    load VCFs and save as hail mt.  
    Save mt as outdir/gatk_unprocessed.mt

    ### Config fields
    ```
    step1.gatk_vcf_header_infile
    step1.gatk_vcf_indir
    step1.gatk_mt_outfile
    ```
    '''
    indir, header, outfile = (
        config['step1']['gatk_vcf_indir'], 
        config['step1'].get('gatk_vcf_header_infile'), # optional
        config['step1']['gatk_mt_outfile']
    )

    objects = hl.utils.hadoop_ls(path_spark(indir))
    
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
        
    mt_out_file = path_spark(outfile)
    print(f"Saving as hail mt to {mt_out_file}")
    mt.write(mt_out_file, overwrite=True)


def main():
    #set up input variables
    config = parse_config()
    
    #initialise hail
    tmp_dir = config['general']['tmp_dir']
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    #load VCFs
    load_vcfs_to_mt(config)


if __name__ == '__main__':
    main() 
