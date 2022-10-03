#export to VCF after filtering 
#TODO fix header. decide what fields are needed in the VCF
import hail as hl
import pyspark
from wes_qc.utils.utils import parse_config


def export_vcfs(mtfile: str, vcf_dir: str):
    '''
    Export VCFs after filtering
    :param str mtfile: Input MatrixTable file
    :param str vcf_dir: Directory for output VCFs
    '''
    metadata = {
        'info': {'consequence': {'Description': 'Most severe consequence from VEP104',
                   'Number': 'A',
                   'Type': 'String'},
                   'rf_score': {'Description': 'Score from variant QC randon forest',
                   'Number': 'A',
                   'Type': 'Float'},
                   'rf_bin': {'Description': 'Bin from variant QC random forest',
                   'Number': 'A',
                   'Type': 'Integer'}
                }
    }

    mt = hl.read_matrix_table(mtfile)

    chroms=[*range(1,23),"X","Y"]
    chromosomes=["chr"+ str(chr) for chr in chroms]
    for chromosome in chromosomes:
        print("Exporting " + chromosome)
        mt_chrom=mt.filter_rows(mt.locus.contig==chromosome)
        outfile = vcf_dir + chromosome + "_filtered.vcf.bgz"
        hl.export_vcf(mt_chrom, outfile, metadata = metadata)


def main():
    #set up
    inputs = parse_config()
    mtdir = inputs['matrixtables_lustre_dir']
    filtered_vcf_dir = inputs['vcf_output_dir']

    # initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    mtfile_filtered = mtdir + "mt_after_var_qc_hard_filter_gt.mt"
    export_vcfs(mtfile_filtered, filtered_vcf_dir)


if __name__ == '__main__':
    main()