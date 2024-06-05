# Load GATK VCFs into hail and save as matrixtable
import hail as hl
import pyspark
import yaml
import os
import sys
import argparse
from utils.utils import parse_config

def get_options():
    '''
    Get options from the command line
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--filter",
        help="filter control samples", action="store_true")
    args = parser.parse_args()

    return args

def load_vcfs_to_mt(indir, outdir, tmp_dir, header, control_list, con_filter):
    '''
    load VCFs and save as hail mt
    '''
    objects = hl.utils.hadoop_ls(indir)
    vcfs = [vcf["path"] for vcf in objects if (vcf["path"].startswith("file") and vcf["path"].endswith("vcf.gz"))]
    print("Loading VCFs")
    #create and save MT
    mt = hl.import_vcf(vcfs, array_elements_required=False, force_bgz=True, header_file=None)
    #separate mt
    if con_filter:
        with open(control_list, 'r') as file:
            control_set = set([line.strip() for line in file])
        con_mt = mt.filter_cols(hl.literal(control_set).contains(mt.s))
        mt = mt.filter_cols(~hl.literal(control_set).contains(mt.s))
        mt_con_out_file = outdir + "gatk_unprocessed.control_samples.mt"
        con_mt.write(mt_con_out_file, overwrite=True)
    print("Saving as hail mt")
    mt_out_file = outdir + "gatk_unprocessed.mt"
    mt.write(mt_out_file, overwrite=True)


def main():
    #set up input variables
    args = get_options()
    inputs = parse_config()
    vcf_header = inputs['gatk_vcf_header']
    import_vcf_dir = inputs['gatk_import_lustre_dir']
    mtdir = inputs['matrixtables_lustre_dir']
    control_list = inputs['control_samples']
    #initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    #load VCFs
    load_vcfs_to_mt(import_vcf_dir, mtdir, tmp_dir, vcf_header, control_list, args.filter)


if __name__ == '__main__':
    main() 