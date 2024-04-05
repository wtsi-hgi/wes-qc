# export to VCF after annotating with a range of hard filters
import os.path
import hail as hl
import pyspark
import argparse
from utils.utils import parse_config

fail_string = 'FAIL'
outside_bait_string = 'OUTSIDE_BAIT'
pass_medium_string = 'PASS_MEDIUM'
pass_stringent_string = 'PASS_STRINGENT'


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
    mt = mt.filter_rows(mt.info.fraction_pass_relaxed_filters > 0)

    # drop unwanted fields
    mt = mt.annotate_rows(info=mt.info.drop('AC', 'AN', 'AF'))
    for filter_level in  ('relaxed', 'medium', 'stringent'):
        for metric in ('AN', 'AC', 'AC_Hom', 'AC_Het'):
            filter_metric = f'{filter_level}_{metric}'
            mt = mt.annotate_rows(info=mt.info.annotate(
                **{filter_metric: mt[filter_metric]}
            ))
            mt = mt.drop(mt[filter_metric])

    mt = mt.drop(
        mt.a_index, mt.was_split, mt.variant_qc,
        mt.stringent_pass_count, mt.medium_pass_count, mt.relaxed_pass_count,
        mt.adj, mt.assigned_pop, mt.sum_AD
    )

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
        str(hard_filters['missingness'])

    medium_filters = "SNPs: RF bin<=" + \
        str(hard_filters['snp']['medium']['bin']) + " & DP>=" + \
        str(hard_filters['snp']['medium']['dp']) + " & GQ>=" + \
        str(hard_filters['snp']['medium']['gq']) + " & HetAB>=" + \
        str(hard_filters['snp']['medium']['ab']) + ", Indels: RF bin<=" + \
        str(hard_filters['indel']['medium']['bin']) + " & DP>=" + \
        str(hard_filters['indel']['medium']['dp']) + " & GQ>=" + \
        str(hard_filters['indel']['medium']['gq']) + " & HetAB>=" + \
        str(hard_filters['indel']['medium']['ab']) + " & Call Rate>=" + \
        str(hard_filters['missingness'])

    relaxed_filters = "SNPs: RF bin<=" + \
        str(hard_filters['snp']['relaxed']['bin']) + " & DP>=" + \
        str(hard_filters['snp']['relaxed']['dp']) + " & GQ>=" + \
        str(hard_filters['snp']['relaxed']['gq']) + " & HetAB>=" + \
        str(hard_filters['snp']['relaxed']['ab']) + ", Indels: RF bin<=" + \
        str(hard_filters['indel']['relaxed']['bin']) + " & DP>=" + \
        str(hard_filters['indel']['relaxed']['dp']) + " & GQ>=" + \
        str(hard_filters['indel']['relaxed']['gq']) + " & HetAB>=" + \
        str(hard_filters['indel']['relaxed']['ab']) + " & Call Rate>=" + \
        str(hard_filters['missingness'])

    metadata = {
        'info': {
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
                                 'Type': 'Integer'},
            'medium_AN': {'Description': 'Total number of alleles in called genotypes',
                          'Number': '1',
                          'Type': 'Integer'},
            'medium_AC': {'Description': 'Allele count in genotypes',
                          'Number': 'A',
                          'Type': 'Integer'},
            'medium_AC_Hom': {'Description': 'Allele counts in homozygous genotypes',
                              'Number': 'A',
                              'Type': 'Integer'},
            'medium_AC_Het': {'Description': 'Allele counts in heterozygous genotypes',
                              'Number': 'A',
                              'Type': 'Integer'},
            'fraction_pass_stringent_filters': {'Description': 'Fraction of genotypes which pass stringent hard filters ' + stringent_filters,
                                                'Number': 'A',
                                                'Type': 'Float'},
            'fraction_pass_medium_filters': {'Description': 'Fraction of genotypes which pass medium hard filters ' + medium_filters,
                                             'Number': 'A',
                                             'Type': 'Float'},
            'fraction_pass_relaxed_filters': {'Description': 'Fraction of genotypes which pass relaxed hard filters ' + relaxed_filters,
                                              'Number': 'A',
                                              'Type': 'Float'},
            'rf_bin': {'Description': 'Variant QC random forest bin, model id ' + run_hash,
                       'Number': 'A',
                       'Type': 'Integer'},
            'rf_score': {'Description': 'Variant QC random forest score, model id ' + run_hash,
                         'Number': 'A',
                         'Type': 'Float'},
            'CSQ': {'Description': 'Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|UNIPROT_ISOFORM|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|gnomADe_AF|gnomADe_AFR_AF|gnomADe_AMR_AF|gnomADe_ASJ_AF|gnomADe_EAS_AF|gnomADe_FIN_AF|gnomADe_NFE_AF|gnomADe_OTH_AF|gnomADe_SAS_AF|gnomADg_AF|gnomADg_AFR_AF|gnomADg_AMI_AF|gnomADg_AMR_AF|gnomADg_ASJ_AF|gnomADg_EAS_AF|gnomADg_FIN_AF|gnomADg_MID_AF|gnomADg_NFE_AF|gnomADg_OTH_AF|gnomADg_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|PHENOTYPES|Conservation',
                    'Number': 'A',
                    'Type': 'String'},
            'consequence': {'Description': 'Most severe consequence from VEP110',
                            'Number': 'A',
                            'Type': 'String'},
            'gene': {'Description': 'Gene affected by the most severe consequence from VEP110',
                     'Number': 'A',
                     'Type': 'String'},
            'hgnc_id': {'Description': 'HGNC id of the gene affected by the most severe consequence from VEP110',
                        'Number': 'A',
                        'Type': 'String'}
        },
        'format': {'HetAB': {'Description': 'Hetrozygous allele balance',
                             'Number': 'A',
                             'Type': 'Float'},
                   'stringent_filters': {'Description': 'Pass/fail stringent hard filters ' + stringent_filters,
                                         'Number': 'A',
                                         'Type': 'String'},
                   'medium_filters': {'Description': 'Pass/fail hard medium filters ' + medium_filters,
                                      'Number': 'A',
                                      'Type': 'String'},
                   'relaxed_filters': {'Description': 'Pass/fail relaxed hard filters ' + relaxed_filters,
                                       'Number': 'A',
                                       'Type': 'String'}

                   },
    }

    #export per chromosome
    chroms = [*range(1, 23), "X", "Y"]
    chromosomes=["chr"+ str(chr) for chr in chroms]
    for chromosome in chromosomes:
        print("Exporting " + chromosome)
        mt_chrom = mt.annotate_globals(chromosome=chromosome)
        mt_chrom = mt_chrom.filter_rows(mt_chrom.locus.contig == mt_chrom.chromosome)
        outfile = os.path.join(filtered_vcf_dir,  f"{chromosome}_hard_filters.vcf.bgz")
        hl.export_vcf(mt_chrom, outfile, metadata=metadata)


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
