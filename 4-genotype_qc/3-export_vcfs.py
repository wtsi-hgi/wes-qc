#export to VCF after filtering 
#TODO fix header. decide what fields are needed in the VCF
import os
import hail as hl
import pyspark
from tempfile import NamedTemporaryFile
from utils.utils import parse_config, path_spark, path_local


def export_vcfs(mtfile: str, vcf_dir: str, config: dict):
    '''
    Export VCFs after filtering
    :param str mtfile: Input MatrixTable file
    :param str vcf_dir: Directory for output VCFs
    '''
    metadata = {
        'info': {
            'consequence': {
                'Description': 'Most severe consequence from VEP110',
                'Number': 'A',
                'Type': 'String'
            },
            'rf_score': {
                'Description': 'Score from variant QC randon forest',
                'Number': 'A',
                'Type': 'Float'
            },
            'rf_bin': {
                'Description': 'Bin from variant QC random forest',
                'Number': 'A',
                'Type': 'Integer'
            }
        }
    }

    mt = hl.read_matrix_table(mtfile)
    #drop boolean format field adj
    mt = mt.drop(mt.adj, mt.hard_filters)

    header_lines = [
        '##source=SLEmap_2024',
        '##genotypeFilters=<bin=86;DP=5;GQ=20;HetAB=0.2>'
    ]

    header = NamedTemporaryFile(dir=path_local(config['tmp_dir']), delete=False)
    with open(header.name, 'w') as header_file:
        header_file.write("\n".join(header_lines))

    chroms = [*range(1, 23), "X", "Y"]
    chromosomes = ["chr" + str(chr) for chr in chroms]
    for chromosome in chromosomes:
        print("Exporting " + chromosome)
        mt_chrom=mt.filter_rows(mt.locus.contig == chromosome)
        outfile = os.path.join(vcf_dir, chromosome + "_filtered.vcf.bgz")
        hl.export_vcf(mt_chrom, outfile, metadata=metadata, append_to_header=path_spark(header.name))

    os.remove(header.name)


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

    mtfile_filtered = os.path.join(mtdir, "mt_after_var_qc_hard_filter_gt.mt")
    export_vcfs(mtfile_filtered, filtered_vcf_dir, config=inputs)


if __name__ == '__main__':
    main()
