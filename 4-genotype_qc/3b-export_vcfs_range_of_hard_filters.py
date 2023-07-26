# export to VCF after annotating with a range of hard filters
import hail as hl
import pyspark
import argparse
from utils.utils import parse_config


def get_options():
    '''
    Get options from the command line
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("--runhash", help="RF run hash")
    args = parser.parse_args()
    if not args.runhash:
        print("--runhash must be specified")
        exit(1)

    return args


def export_vcfs(mtfile: str, filtered_vcf_dir: str, hard_filters: dict, run_hash: str):
    '''
    Export VCFs annotated with a range of hard filters
    :param str mtfile: matrixtable file
    :param str filtered_vcf_dir: output directory for VCFs
    :param dict hard_filters: details of all sets of filters
    :param str run_hash: random forest run hash used
    '''
    mt = hl.read_matrix_table(mtfile)
    # filter to remove rows where all variants fail the most relaxed filters
    mt = mt.filter_rows(mt.info.fraction_pass_stringent_filters > 0)
    mt = mt.filter_entries(mt.stringent_filters == 'Pass')
    # drop unwanted fields
    mt = mt.drop(mt.a_index, mt.was_split, mt.stringent_pass_count, mt.medium_pass_count,
                 mt.adj, mt.assigned_pop, mt.sum_AD,
                 mt.medium_filters, mt.stringent_filters,
                 mt.medium_AN, mt.medium_AC, mt.medium_AC_Hom, mt.medium_AC_Het,
                 mt.relaxed_AN, mt.relaxed_AC, mt.relaxed_AC_Hom, mt.relaxed_AC_Het)
    mt = mt.annotate_rows(info=mt.info.drop('CSQ', 'consequence', 'gene', 'hgnc_id',
                                            'fraction_pass_medium_filters', 'fraction_pass_relaxed_filters'))
    mt = mt.annotate_rows(info=mt.info.annotate(stringent_AN=mt.stringent_AN,
                                                stringent_AC=mt.stringent_AC,
                                                stringent_AC_Hom=mt.stringent_AC_Hom,
                                                stringent_AC_Het=mt.stringent_AC_Het))
    # info for header
    stringent_filters = "SNPs: RF bin<=" + \
        str(hard_filters['snp']['stringent']['bin']) + " & DP>=" + \
        str(hard_filters['snp']['stringent']['dp']) + " & GQ>=" + \
        str(hard_filters['snp']['stringent']['gq']) + " & HetAB>=" + \
        str(hard_filters['snp']['stringent']['ab']) + ", Indels: RF bin<=" + \
        str(hard_filters['indel']['stringent']['bin']) + " & DP>=" + \
        str(hard_filters['indel']['stringent']['dp']) + " & GQ>=" + \
        str(hard_filters['indel']['stringent']['gq']) + " & HetAB>=" + \
        str(hard_filters['indel']['stringent']['ab']) + " & Call Rate>=" + \
        str(hard_filters['missingness']['stringent'])

    medium_filters = "SNPs: RF bin<=" + \
        str(hard_filters['snp']['medium']['bin']) + " & DP>=" + \
        str(hard_filters['snp']['medium']['dp']) + " & GQ>=" + \
        str(hard_filters['snp']['medium']['gq']) + " & HetAB>=" + \
        str(hard_filters['snp']['medium']['ab']) + ", Indels: RF bin<=" + \
        str(hard_filters['indel']['medium']['bin']) + " & DP>=" + \
        str(hard_filters['indel']['medium']['dp']) + " & GQ>=" + \
        str(hard_filters['indel']['medium']['gq']) + " & HetAB>=" + \
        str(hard_filters['indel']['medium']['ab'])     

    relaxed_filters = "SNPs: RF bin<=" + \
        str(hard_filters['snp']['relaxed']['bin']) + " & DP>=" + \
        str(hard_filters['snp']['relaxed']['dp']) + " & GQ>=" + \
        str(hard_filters['snp']['relaxed']['gq']) + " & HetAB>=" + \
        str(hard_filters['snp']['relaxed']['ab']) + ", Indels: RF bin<=" + \
        str(hard_filters['indel']['relaxed']['bin']) + " & DP>=" + \
        str(hard_filters['indel']['relaxed']['dp']) + " & GQ>=" + \
        str(hard_filters['indel']['relaxed']['gq']) + " & HetAB>=" + \
        str(hard_filters['indel']['relaxed']['ab']) 

    metadata = {
        'filter': {
            'ExcessHet': {'Description': 'ExcessHet > 54.69'},
            'LowQual': {'Description': 'Low quality'}
        },
        'format': {
            'HetAB': {'Description': 'Heterozygous allele balance',
                      'Number': 'A',
                      'Type': 'Float'},
        },
        'info': {
            'fraction_pass_stringent_filters': {'Description': 'Fraction of genotypes which pass stringent hard filters ' + stringent_filters,
                                                'Number': 'A',
                                                'Type': 'Float'},
            'rf_score': {'Description': 'Variant QC random forest score, model id ' + run_hash,
                         'Number': 'A',
                         'Type': 'Float'},
            'rf_bin': {'Description': 'Variant QC random forest bin, model id ' + run_hash,
                       'Number': 'A',
                       'Type': 'Integer'},
            'stringent_AN': {'Description': 'Total number of alleles in called genotypes',
                             'Number': '1',
                             'Type': 'Integer'},
            'stringent_AC': {'Description': 'Allele count in genotypes',
                             'Number': 'A',
                             'Type': 'Integer'},
            'stringent_AC_Hom': {'Description': 'Allele counts in homozygous genotypes',
                                 'Number': 'A',
                                 'Type': 'Integer'},
            'stringent_AC_Het': {'Description': 'Allele counts in heterozygous genotypes',
                                 'Number': 'A',
                                 'Type': 'Integer'}
        }
    }

    #export per chromosome
    # chroms = [*range(1,23),"X","Y"]
    chroms = [22]
    chromosomes = ["chr"+ str(chr) for chr in chroms]
    for chromosome in chromosomes:
        print("Exporting " + chromosome)
        mt_chrom = mt.annotate_globals(chromosome=chromosome)
        mt_chrom = mt_chrom.filter_rows(mt_chrom.locus.contig == mt_chrom.chromosome)
        outfile = filtered_vcf_dir + chromosome + "_hard_filters.vcf.bgz"
        hl.export_vcf(mt_chrom, outfile, metadata = metadata)


def main():
    # set up
    inputs = parse_config()
    args = get_options()
    mtdir = inputs['matrixtables_lustre_dir']
    hard_filters = inputs['hard_filters']
    filtered_vcf_dir = inputs['vcf_output_dir']

    # initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    mtfile = mtdir + "mt_hard_filter_combinations.mt"
    export_vcfs(mtfile, filtered_vcf_dir, hard_filters, args.runhash)


if __name__ == '__main__':
    main()
