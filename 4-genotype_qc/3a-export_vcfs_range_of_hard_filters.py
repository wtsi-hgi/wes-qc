# export to VCF after annotating with a range of hard filters
import hail as hl
import pyspark
import argparse
from wes_qc.utils.utils import parse_config


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
    mt = mt.drop(mt.a_index, mt.was_split, mt.stringent_pass_count, mt.stringent_fail_count,
                 mt.medium_pass_count, mt.medium_fail_count, mt.relaxed_pass_count, mt.relaxed_fail_count,
                 mt.batch, mt.sequencing_location, mt.adj, mt.assigned_pop, mt.sum_AD)
    # info for header
    stringent_filters = "SNPs: RF bin<=" + \
        str(hard_filters['snp']['stringent']['bin']) + " & DP>=" + \
        str(hard_filters['snp']['stringent']['dp']) + " & GQ>=" + \
        str(hard_filters['snp']['stringent']['gq']) + " & HetAB>=" + \
        str(hard_filters['snp']['stringent']['ab']) + ", Indels: RF bin<=" + \
        str(hard_filters['indel']['stringent']['bin']) + " & DP>=" + \
        str(hard_filters['indel']['stringent']['dp']) + " & GQ>=" + \
        str(hard_filters['indel']['stringent']['gq']) + " & HetAB>=" + \
        str(hard_filters['indel']['stringent']['ab']) 

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
        str(hard_filters['indel']['relaxd']['ab']) 

    metadata = {
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
        'info': {
            'fraction_pass_stringent_filters': {'Description': 'Fraction of genotypes which pass stringent hard filters ' + stringent_filters,
                                                'Number': 'A',
                                                'Type': 'Float'},
            'fraction_pass_medium_filters': {'Description': 'Fraction of genotypes which pass medium hard filters ' + medium_filters,
                                             'Number': 'A',
                                             'Type': 'Float'},
            'fraction_pass_relaxed_filters': {'Description': 'Fraction of genotypes which pass relaxed hard filters ' + relaxed_filters,
                                              'Number': 'A',
                                              'Type': 'Float'},
            'rf_score': {'Description': 'Variant QC random forest score, model id ' + run_hash,
                         'Number': 'A',
                         'Type': 'Float'},
            'rf_bin': {'Description': 'Variant QC random forest bin, model id ' + run_hash,
                       'Number': 'A',
                       'Type': 'Integer'},
            'CSQ': {'Description': 'Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|UNIPROT_ISOFORM|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|SpliceRegion|GeneSplicer|existing_InFrame_oORFs|existing_OutOfFrame_oORFs|existing_uORFs|five_prime_UTR_variant_annotation|five_prime_UTR_variant_consequence|CADD_PHRED|CADD_RAW|Ensembl_transcriptid|LRT_pred|MutationTaster_pred|Polyphen2_HDIV_pred|Polyphen2_HVAR_pred|SIFT_pred|Uniprot_acc|VEP_canonical|DisGeNET_PMID|DisGeNET_SCORE|PHENOTYPES|Conservation|LoF|LoF_filter|LoF_flags|LoF_info|REVEL|SpliceAI_pred_DP_AG|SpliceAI_pred_DP_AL|SpliceAI_pred_DP_DG|SpliceAI_pred_DP_DL|SpliceAI_pred_DS_AG|SpliceAI_pred_DS_AL|SpliceAI_pred_DS_DG|SpliceAI_pred_DS_DL|SpliceAI_pred_SYMBOL',
                    'Number': 'A',
                    'Type': 'String'},
            'consequence': {'Description': 'Most severe consequence from VEP104',
                            'Number': 'A',
                            'Type': 'String'},
            'gene': {'Description': 'Gene affected by the most severe consequence from VEP104',
                     'Number': 'A',
                     'Type': 'String'},
            'hgnc_id': {'Description': 'HGNC id of the gene affected by the most severe consequence from VEP104',
                        'Number': 'A',
                        'Type': 'String'}
        }
    }

    #export per chromosome
    chroms=[*range(1,23),"X","Y"]
    chromosomes=["chr"+ str(chr) for chr in chroms]
    for chromosome in chromosomes:
        print("Exporting " + chromosome)
        mt_chrom=mt.filter_rows(mt.locus.contig==chromosome)
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
