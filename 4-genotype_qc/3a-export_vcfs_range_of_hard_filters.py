# export to VCF after annotating with a range of hard filters
import os.path
import hail as hl
import pyspark
import argparse
from pathlib import Path
from typing import Union
from utils.utils import parse_config
from utils.utils import remove_samples
from utils.utils import rm_mt

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


def annotate_variants(mt: hl.MatrixTable, mt_filtered: hl.MatrixTable) -> hl.MatrixTable:
    mt_filtered_locus = mt_filtered.key_rows_by(mt_filtered.locus)
    mt = mt.annotate_rows(
        filters=hl.case()
                  .when(~hl.is_defined(mt_filtered_locus.index_rows(mt.locus)), {outside_bait_string})
                  .when(mt_filtered.index_rows(mt.row_key).info.fraction_pass_stringent_filters > 0, {pass_stringent_string})
                  .when(mt_filtered.index_rows(mt.row_key).info.fraction_pass_medium_filters > 0, {pass_medium_string})
                  .default({fail_string}),
        info=mt.info.annotate(rf_bin=mt_filtered.index_rows(mt.row_key).info.rf_bin),
    )
    return mt


def annotate_genotypes(mt: hl.MatrixTable, hard_filters) -> hl.MatrixTable:
    mt = mt.annotate_entries(
        HetAB=hl.case()
            .when(mt.sum_AD == 0, hl.missing('float64'))
            .when(mt.GT.is_het(), hl.min(mt.AD[1]) / mt.sum_AD)
            .or_missing()
    )

    stringent_condition = (
            (
                    (hl.is_snp(mt.alleles[0], mt.alleles[1]))
                    & (mt.info.rf_bin <= hard_filters['snp']['stringent']['bin'])
                    & (mt.DP >= hard_filters['snp']['stringent']['dp'])
                    & (mt.GQ >= hard_filters['snp']['stringent']['gq'])
                    & (
                            (mt.GT.is_het() & (mt.HetAB >= hard_filters['snp']['stringent']['ab'])) |
                            (mt.GT.is_hom_ref()) |
                            (mt.GT.is_hom_var())
                    )

            ) |
            (
                    (hl.is_indel(mt.alleles[0], mt.alleles[1]))
                    & (mt.info.rf_bin <= hard_filters['indel']['stringent']['bin'])
                    & (mt.DP >= hard_filters['indel']['stringent']['dp'])
                    & (mt.GQ >= hard_filters['indel']['stringent']['gq'])
                    & (
                            (mt.GT.is_het() & (mt.HetAB >= hard_filters['indel']['stringent']['ab'])) |
                            (mt.GT.is_hom_ref()) |
                            (mt.GT.is_hom_var())
                    )

            )
    )

    medium_condition = (
            (
                    (hl.is_snp(mt.alleles[0], mt.alleles[1]))
                    & (mt.info.rf_bin <= hard_filters['snp']['medium']['bin'])
                    & (mt.DP >= hard_filters['snp']['medium']['dp'])
                    & (mt.GQ >= hard_filters['snp']['medium']['gq'])
                    & (
                            (mt.GT.is_het() & (mt.HetAB >= hard_filters['snp']['medium']['ab'])) |
                            (mt.GT.is_hom_ref()) |
                            (mt.GT.is_hom_var())
                    )

            ) |
            (
                    (hl.is_indel(mt.alleles[0], mt.alleles[1]))
                    & (mt.info.rf_bin <= hard_filters['indel']['medium']['bin'])
                    & (mt.DP >= hard_filters['indel']['medium']['dp'])
                    & (mt.GQ >= hard_filters['indel']['medium']['gq'])
                    & (
                            (mt.GT.is_het() & (mt.HetAB >= hard_filters['indel']['medium']['ab'])) |
                            (mt.GT.is_hom_ref()) |
                            (mt.GT.is_hom_var())
                    )

            )
    )

    mt = mt.annotate_entries(
        status=hl.case()
        .when(mt.filters.contains(outside_bait_string), outside_bait_string)
        .when(stringent_condition, pass_stringent_string)
        .when(medium_condition, pass_medium_string)
        .default(fail_string)
    )
    return mt


def annotate_ac(mt: hl.MatrixTable, filter_name: str) -> hl.MatrixTable:
    if filter_name == 'stringent':
        mt_filtered = mt.filter_entries(mt.status == pass_stringent_string)
    elif filter_name == 'medium':
        mt_filtered = mt.filter_entries((mt.status == pass_stringent_string) | (mt.status == pass_medium_string))
    else:
        raise ValueError('filter_name must be `stringent` or `medium`')

    mt_filtered = hl.variant_qc(mt_filtered, name='vqc')

    vqc = mt_filtered.index_rows(mt.row_key).vqc
    annotation = {
        f'{filter_name}_AN': vqc.AN,
        f'{filter_name}_AC': vqc.AC[1:],
        f'{filter_name}_AC_Hom': 2 * vqc.homozygote_count[1:],
        f'{filter_name}_AC_Het': vqc.AC[1:] - 2*vqc.homozygote_count[1:]
    }
    mt = mt.annotate_rows(info=mt.info.annotate(**annotation))
    return mt


def export_vcfs(mtfile: str, mtrawfile: str, filtered_vcf_dir: str, hard_filters: dict, run_hash: str,
                samples_failing_qc_path: Union[Path, str], samples_to_remove_path: Union[Path, str], chroms):
    '''
    Export VCFs annotated with a range of hard filters
    :param str mtfile: matrixtable file
    :param str filtered_vcf_dir: output directory for VCFs
    :param dict hard_filters: details of all sets of filters
    :param str run_hash: random forest run hash used
    '''
    mt_qc = hl.read_matrix_table(mtfile)
    mt = hl.read_matrix_table(mtrawfile)
    mt = mt.annotate_entries(sum_AD=hl.sum(mt.AD))
    mt = hl.split_multi_hts(mt)
    mt = mt.annotate_rows(info=mt.info.drop('AC', 'AC_adj', 'AN'))

    # exclude samples failing QC
    mt = remove_samples(mt, samples_failing_qc_path)

    # exclude extra samples
    mt = remove_samples(mt, samples_to_remove_path, sample_column_name='ID')

    # flag variants
    mt = annotate_variants(mt, mt_qc)

    # flag genotypes
    mt = annotate_genotypes(mt, hard_filters)
    mt = mt.drop(mt.sum_AD)

    # add extra INFO annotation
    for filter_name in ('stringent', 'medium'):
        mt = annotate_ac(mt, filter_name)

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
        str(hard_filters['indel']['medium']['ab']) + " & Call Rate>=" + \
        str(hard_filters['missingness']['medium'])

    metadata = {
        'info': {
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
                              'Type': 'Integer'}
                 },
        'format': {'HetAB': {'Description': 'Heterozygous allele balance',
                             'Number': 'A',
                             'Type': 'Float'},
                   'status': {'Description': 'One of PASS_MEDIUM or PASS_STRINGENT or FAIL or OUTSIDE_BAIT',
                              'Number': 'A',
                              'Type': 'String'}
                   },
        'filter': {
            outside_bait_string: {'Description': 'Variant is outside baits +-50bp'},
            fail_string: {'Description': 'Variant fails all filters'},
            pass_medium_string: {
                'Description': f"Variant passes medium filters: RF bin<={hard_filters['snp']['medium']['bin']} & Call Rate>={hard_filters['missingness']['medium']}"
            },
            pass_stringent_string: {
                'Description': f"Variant passes stringent filters: RF bin<={hard_filters['snp']['stringent']['bin']} & Call Rate>={hard_filters['missingness']['stringent']}"
            }
        }
    }

    temppath = os.path.join(filtered_vcf_dir, 'tmp.mt')
    mt = mt.checkpoint(temppath)

    #export per chromosome
    #chroms=[*range(1,23),"X","Y"]
    #chroms = ['X']
    chromosomes=["chr"+ str(chr) for chr in chroms]
    for chromosome in chromosomes:
        print("Exporting " + chromosome)
        mt_chrom = mt.annotate_globals(chromosome=chromosome)
        mt_chrom = mt_chrom.filter_rows(mt_chrom.locus.contig == mt_chrom.chromosome)
        outfile = os.path.join(filtered_vcf_dir,  f"{chromosome}_hard_filters.vcf.bgz")
        hl.export_vcf(mt_chrom, outfile, metadata=metadata)

    rm_mt(temppath)


def main():
    # set up
    inputs = parse_config()
    args = get_options()
    mtdir = inputs['matrixtables_lustre_dir']
    hard_filters = inputs['hard_filters']
    filtered_vcf_dir = inputs['vcfraw_output_dir']
    annotation_dir = Path(inputs['annotation_lustre_dir'])

    # initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    mtraw_file = mtdir + 'gatk_unprocessed.mt'
    mtfile = mtdir + "mt_hard_filter_combinations.mt"
#    mtfile = mtdir + "mt_hard_filter_combinations.chrX.mt"
    samples_failing_qc = annotation_dir / 'samples_failing_qc.tsv'
    samples_to_remove = inputs['exclude_samples']
    filtered_vcf_dir = Path(filtered_vcf_dir).parent / 'annotated_vcfs'
    chroms=[*range(1,23),"Y"]
    export_vcfs(mtfile, mtraw_file, str(filtered_vcf_dir), hard_filters, args.runhash,
                samples_failing_qc, samples_to_remove, chroms)
    mtfile = mtdir + "mt_hard_filter_combinations.chrX.mt"
    chroms=["X"]
    export_vcfs(mtfile, mtraw_file, str(filtered_vcf_dir), hard_filters, args.runhash,
                samples_failing_qc, samples_to_remove, chroms)

if __name__ == '__main__':
    main()
