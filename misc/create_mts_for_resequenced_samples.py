#create mt for the 112 resequenced samples to feed into variant QC
#also create a MT for the original run of sequencing from these to use in variant QC as a comparison
import hail as hl
import pyspark
from wes_qc.utils.utils import parse_config


def import_resequenced_vcfs(import_resequenced_vcf_dir: str, vcf_header: str, resequenced_mt_file: str):
    '''
    :param str import_resequenced_vcf_dir: directiry containing VCFs to import
    :param str vcf_header: path to VCF header file
    :param str resequenced_mt_file: output file name
    '''
    objects = hl.utils.hadoop_ls(import_resequenced_vcf_dir)
    vcfs = [vcf["path"] for vcf in objects if (vcf["path"].startswith("file") and vcf["path"].endswith("vcf.gz"))]
    print("Loading resequenced VCFs")
    mt = hl.import_vcf(vcfs, array_elements_required=False, force_bgz=True, header_file = vcf_header)
    print("Saving as hail mt")
    mt.write(resequenced_mt_file, overwrite=True)


def subset_mt_by_sample(original_mtfile: str, samplefile: str, output_mtfile: str):
    '''
    :param str original_mtfile: original mt file to subset by sample
    :param str samplefile: samples to include
    :param str output_mtfile: file for output matrixtable
    '''
    mt = hl.read_matrix_table(original_mtfile)
    sample_ht = hl.import_table(samplefile, delimiter="\t").key_by('s')
    mt = mt.filter_cols(hl.is_defined(sample_ht[mt.s]))
    mt.write(output_mtfile, overwrite = True)


def main():
    #set up input variables
    inputs = parse_config()
    # mtdir = inputs['load_matrixtables_lustre_dir']
    annot_dir = inputs['annotation_lustre_dir']
    mtdir = inputs['matrixtables_lustre_dir']

    #initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    #import new sequencing data for the 112 resequenced samples
    import_resequenced_vcf_dir = 'file:///lustre/scratch123/qc/repeat_samples_gatk/'
    vcf_header = import_resequenced_vcf_dir + "header.txt"
    resequenced_mt_file = mtdir + 'resequenced_samples_gatk_unprocessed.mt'
    import_resequenced_vcfs(import_resequenced_vcf_dir, vcf_header, resequenced_mt_file)

    #subset original mt to just the resequenced samples
    original_mtfile = mtdir + "gatk_unprocessed.mt"
    samplefile = import_resequenced_vcf_dir + "resequenced_samples.txt"
    high_snp_samples_mtfile = mtdir + "gatk_unprocessed_high_snp_samples.mt"
    subset_mt_by_sample(original_mtfile, samplefile, high_snp_samples_mtfile)

if __name__ == '__main__':
    main()
