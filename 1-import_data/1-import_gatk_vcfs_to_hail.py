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
    
    vcf_names = r'\file://.*vcf.b?gz\'
    
    # TODO: add .bgz support
    vcfs = [vcf["path"] for vcf in objects if vcf_names.match(vcf["path"])]
    #create and save MT

    # TODO: make header file parameter optional
    if header:
        print("Loading VCFs with header")
        mt = hl.import_vcf(vcfs, array_elements_required=False, force_bgz=True, header_file = header)
    else:
        print("Loading VCFs WITHOUT header")
        mt = hl.import_vcf(vcfs, array_elements_required=False, force_bgz=True)
        
    print("Saving as hail mt")
    # Using prefix mk43 to discern between my and Dmytro's data, as we are sharing an mt directory
    mt_out_file = os.path.join(outdir, "mk43_gatk_unprocessed.mt")
    mt.write(mt_out_file, overwrite=True)


def main():
    #set up input variables
    inputs = parse_config()
    vcf_header = inputs['gatk_vcf_header']
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
