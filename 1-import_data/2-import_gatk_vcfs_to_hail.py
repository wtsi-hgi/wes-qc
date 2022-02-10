# Load GATK VCFs into hail and save as matrixtable
import hail as hl
import pyspark
import yaml
import os
import sys
from wes_qc.utils.utils import parse_config


def load_vcfs_to_mt(indir, outdir, tmp_dir, header):
    '''
    load VCFs and save as hail mt and as a sparse_mt
    '''
    objects = hl.utils.hadoop_ls(indir)
    vcfs = [vcf["path"] for vcf in objects if (vcf["path"].startswith("file") and vcf["path"].endswith("vcf.gz"))]
    print("Loading VCFs")
    #create and save normal MT
    # mt = hl.import_vcf(vcfs, array_elements_required=False, force_bgz=True, header_file = header)
    # print("Saving as hail mt")
    # mt_out_file = outdir + "gatk_unprocessed.mt"
    # mt.write(mt_out_file, overwrite=True)
    #create sparse mt
    vcfs = ['file:///lustre/scratch123/qc/gatk_vcfs/chr17_64212474_chr17_64531131.gatk.vcf.gz']
    sparse_mt_out_file = outdir + "gatk_unprocessed_sparse.mt"
    hl.experimental.run_combiner(vcfs, out_file=sparse_mt_out_file, tmp_path=tmp_dir, use_exome_default_intervals=True, reference_genome='GRCh38')


def main():
    #set up input variables
    inputs = parse_config()
    vcf_header = inputs['gatk_vcf_header']
    import_vcf_dir = inputs['gatk_import_lustre_dir']
    inputs = parse_config()
    mtdir = inputs['matrixtables_lustre_dir']

    #initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")
    #hl.init(sc=sc, tmp_dir=lustre_dir, local_tmpdir=lustre_dir, default_reference="GRCh38")

    #load VCFs
    load_vcfs_to_mt(import_vcf_dir, mtdir, tmp_dir, vcf_header)


if __name__ == '__main__':
    main() 