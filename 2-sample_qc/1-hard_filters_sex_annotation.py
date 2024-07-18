#apply gnomad's hard filters and impute sex
#input gatk_unprocessed.mt from step 1.1
import hail as hl
import pyspark
from utils.utils import parse_config

import os

def apply_hard_filters(mt: hl.MatrixTable, mtdir: str) -> hl.MatrixTable:
    '''
    Applies hard filters and annotates samples in the filtered set with call rate
    :param MatrixTable mt: MT containing samples to be ascertained for sex
    :param str mtdir: directory output matrix tables are written to
    :return: MatrixTable with hard filtering annotation
    :rtype: MatrixTable
    '''
    print("Applying hard filters")
    filtered_mt_file = os.path.join(mtdir, "mt_hard_filters_annotated.mt") # output

    config = parse_config()['step2']['sex_annotation_hard_filters']

    mt = mt.filter_rows((hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1]) &
        (hl.agg.mean(mt.GT.n_alt_alleles()) / 2 > config['n_alt_alleles_threshold']) &
        (hl.agg.fraction(hl.is_defined(mt.GT)) > config['defined_gt_frac_threshold']))
    mt = mt.annotate_cols(callrate=hl.agg.fraction(hl.is_defined(mt.GT)))
    mt.write(filtered_mt_file, overwrite=True)

    return mt


def impute_sex(mt: hl.MatrixTable, mtdir: str, annotdir: str, male_threshold: float = 0.8, female_threshold: float = 0.5, aaf_threshold: float = 0.05) -> hl.MatrixTable:
    '''
    Imputes sex, exports data, and annotates mt with this data
    :param MatrixTable mt: MT containing samples to be ascertained for sex
    :param str mtdir: directory output matrix tables are written to
    :param str annotdir: directory annotation files are written to
    :return: MatrixTable with imputed sex annotations stashed in column annotation 'sex_check'
    :rtype: MatrixTable
    '''
    print("Imputing sex with male_threshold = " + str(male_threshold) + " and female threshold = " + str(female_threshold))

    #filter to X and select unphased diploid genotypes - no need to filter to X as impute_sex takes care of this
    #mt1 = hl.filter_intervals(mt, [hl.parse_locus_interval('chrX')])
    mt1 = hl.split_multi_hts(mt)
    mtx_unphased = mt1.select_entries(GT=hl.unphased_diploid_gt_index_call(mt1.GT.n_alt_alleles()))
    #imput sex on the unphased diploid GTs
    sex_ht = hl.impute_sex(mtx_unphased.GT, aaf_threshold=aaf_threshold, female_threshold=female_threshold, male_threshold=male_threshold)
    #export
    sex_ht.export(os.path.join(annotdir, 'sex_annotated.sex_check.txt.bgz')) # output
    #annotate input (all chroms) mt with imputed sex and write to file
    sex_colnames = ['f_stat', 'is_female']
    sex_ht = sex_ht.select(*sex_colnames)
    mt = mt.annotate_cols(**sex_ht[mt.col_key])
    sex_mt_file = os.path.join(mtdir, "mt_sex_annotated.mt") # output
    print("Writing to " + sex_mt_file)
    mt.write(sex_mt_file, overwrite=True)

    return mt


def identify_inconsistencies(mt: hl.MatrixTable, mtdir: str, 
                             sex_metadata_file: str,
                             conflicting_sex_report_path: str,
                             fstat_outliers_report_path: str,
                             fstat_low: float = 0.2,
                             fstat_high: float = 0.8
                            ):
    '''
    Find samples where annotated sex conflicts with the sex in our metadata
    Find samples where sex is not annotated
    Find samples where f_stat is between fstat_low and fstat_high
    :param MatrixTable mt: MT containing imputed sex in column 'sex_check'
    :param str mtdir: directory output matrix tables are written to
    :param str sex_metadata_file: INPUT, TODO explain metadata structure and constants
    :param str conflicting_sex_report_path: OUTPUT, TODO
    :param str fstat_outliers_report_path: OUTPUT, TODO
    :param float: fstat_low
    :param float: fstat_high
    '''
    print("Annotating samples with inconsistencies:")
    qc_ht = mt.cols()
    #convert is_female boolean to sex
    sex_expr = (hl.case()
            .when(qc_ht.is_female, "female")
            .when(qc_ht.is_female == False, "male")
            .default("undetermined"))
    qc_ht = qc_ht.annotate(sex=sex_expr).key_by('s')

    #annotate with manifest sex - keyed on ega to match identifiers in matrixtable

    metadata_ht = hl.import_table(sex_metadata_file, delimiter="\t").key_by('accession_number')
    #we only want those from the metadata file where sex is known
    metadata_ht = metadata_ht.filter((metadata_ht.gender == 'Male') | (metadata_ht.gender == 'Female'))

    #annotate the sex-predictions with the manifest sex annotation - need to use a join here
    ht_joined = qc_ht.annotate(manifest_sex = metadata_ht[qc_ht.s].gender)

    #identify samples where imputed sex and manifest sex conflict
    conflicting_sex_ht = ht_joined.filter(((ht_joined.sex == 'male') & (ht_joined.manifest_sex == 'Female')) | (
        (ht_joined.sex == 'female') & (ht_joined.manifest_sex == 'Male')))
    conflicting_sex_ht.export(conflicting_sex_report_path) # output

    #identify samples where f stat is between fstat_low and fstat_high
    f_stat_ht = qc_ht.filter( (qc_ht.f_stat > fstat_low) & (qc_ht.f_stat < fstat_high) )
    f_stat_ht.export(fstat_outliers_report_path) # output
    

def main():
    #set up
    inputs = parse_config()
    #importmtdir = inputs['load_matrixtables_lustre_dir']
    mtdir = inputs['matrixtables_lustre_dir']
    annotdir = inputs['annotation_lustre_dir']
    resourcedir = inputs['resource_dir']

    #initialise hail
    tmp_dir = inputs['tmp_dir']
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    mt_in_file = os.path.join(mtdir, "gatk_unprocessed.mt") # input from 1.1
    print("Reading input matrix")
    mt_unfiiltered = hl.read_matrix_table(mt_in_file)

    #apply hard fitlers
    mt_filtered = apply_hard_filters(mt_unfiiltered, mtdir)

    #impute sex
    # TODO: move male_threshold to config?
    male_threshold, female_threshold, aaf_threshold, fstat_low, fstat_high = (
        inputs['step2']['impute_sex']['male_threshold'], 
        inputs['step2']['impute_sex']['female_threshold'],
        inputs['step2']['impute_sex']['aaf_threshold']
    )
    mt_sex = impute_sex(mt_filtered, mtdir, annotdir, male_threshold, female_threshold, aaf_threshold)

    # TODO: where is this function?
    # annotate_ambiguous_sex(mt_sex, mtdir)

    #identify_inconsistencies

    # TODO: automate via namedtuples or kwargs expansion?
    sex_metadata_file, conflicting_sex_report_path, fstat_outliers_report_path, fstat_low, fstat_high = (
        inputs['step2']['sex_inconsistencies']['sex_metadata_file'],
        inputs['step2']['sex_inconsistencies']['conflicting_sex_report_path'],
        inputs['step2']['sex_inconsistencies']['fstat_outliers_report_path'],
        inputs['step2']['sex_inconsistencies']['fstat_low'],
        inputs['step2']['sex_inconsistencies']['fstat_high'],
    )

    # TODO: do we need such a detailed logging, or a single if (... and ... and ...) will suffice?
    shall_find_inconsistencies = True
    if not os.path.exists(sex_metadata_file):
        print("error: identify_inconsistencies: sex_metadata_file missing")
        shall_find_inconsistencies = False
    if not os.path.exists(conflicting_sex_report_path):
        print("error: identify_inconsistencies: conflicting_sex_report_path missing")
        shall_find_inconsistencies = False
    if not os.path.exists(fstat_outliers_report_path):
        print("error: identify_inconsistencies: fstat_outliers_report_path missing")
        shall_find_inconsistencies = False
    if shall_find_inconsistencies:
        identify_inconsistencies(mt_sex, mtdir, sex_metadata_file, conflicting_sex_report_path, 
                                 fstat_outliers_report_path, fstat_low, fstat_high)
    else:
        print("skip identify_inconsistencies because of previous errors")

if __name__ == '__main__':
    main() 

