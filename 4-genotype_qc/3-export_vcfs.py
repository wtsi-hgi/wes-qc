#export to VCF after filtering 
#TODO fix header. decide what fields are needed in the VCF
import hail as hl
import pyspark
from utils.utils import parse_config, path_spark, path_local


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

    mt = hl.read_matrix_table(path_spark(mtfile))
    #drop boolean format field adj
    mt = mt.drop(mt.adj)

    chroms=[*range(1,23), "X", "Y"]
    chromosomes=["chr"+ str(chr) for chr in chroms]
    for chromosome in chromosomes:
        print("Exporting " + chromosome)
        mt_chrom = mt.filter_rows(mt.locus.contig==chromosome)
        outfile = vcf_dir + chromosome + "_filtered.vcf.bgz"
        hl.export_vcf(mt_chrom, outfile, metadata = metadata)


def main():
    #set up
    config = parse_config()
    mtdir = config['general']['matrixtables_dir']

    # initialise hail
    tmp_dir = config['general']['tmp_dir']
    sc = pyspark.SparkContext.getOrCreate()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38", idempotent=True)

    filtered_vcf_dir = config['step4']['export_vcfs']['vcf_output_dir']
    mtfile_filtered = config['step4']['export_vcfs']['mtfile_filtered']
    export_vcfs(mtfile_filtered, filtered_vcf_dir)


if __name__ == '__main__':
    main()