# Load GATK VCFs into hail and save as matrixtable
import hail as hl
import pyspark
import yaml
import os
import sys
from utils.utils import parse_config


def load_vcfs_to_mt(indir, outdir):
    '''
    load VCFs and save as hail mt
    '''
    objects = hl.utils.hadoop_ls(indir)
    vcfs = [vcf["path"] for vcf in objects if (vcf["path"].startswith("file") and vcf["path"].endswith("vcf.gz"))]
    assert len(vcfs) == 1

    print("Loading VCFs")
    mt = hl.import_vcf(vcfs[0], array_elements_required=False, force_bgz=True)

    print("Saving as hail mt")
    mt_out_file = os.path.join(outdir, "gatk_unprocessed.mt")
    mt.write(mt_out_file, overwrite=True)


def main():
    #set up input variables
    inputs = parse_config()
    import_vcf_dir = inputs['gatk_import_lustre_dir']
    mtdir = inputs['matrixtables_lustre_dir']

    #initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    #load VCFs
    load_vcfs_to_mt(import_vcf_dir, mtdir)


if __name__ == '__main__':
    main()
