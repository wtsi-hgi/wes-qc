#apply gnomad's hard filters and impute sex
#input gatk_unprocessed.mt from step 1.1
import hail as hl
import pyspark
from utils.utils import parse_config

import os


# TODO: change manual path joins to os.path.join() to enhance robustness
# note: os.path.join arguments must not start from '/' unless you meant to start from the fs root

def apply_hard_filters(mt: hl.MatrixTable, mtdir: str) -> hl.MatrixTable:
    '''
    Applies hard filters and annotates samples in the filtered set with call rate
    :param MatrixTable mt: MT containing samples to be ascertained for sex
    :param str mtdir: directory output matrix tables are written to
    :return: MatrixTable with hard filtering annotation
    :rtype: MatrixTable
    '''
    print("Applying hard filters")
    # User prefix
    filtered_mt_file = os.path.join(mtdir, "mk43_mt_hard_filters_annotated.mt") # output

    # TODO: move these number to config
    mt = mt.filter_rows((hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1]) &
        (hl.agg.mean(mt.GT.n_alt_alleles()) / 2 > 0.001) &
        (hl.agg.fraction(hl.is_defined(mt.GT)) > 0.99))
    mt = mt.annotate_cols(callrate=hl.agg.fraction(hl.is_defined(mt.GT)))
    mt.write(filtered_mt_file, overwrite=True)

    return mt


def impute_sex(mt: hl.MatrixTable, mtdir: str, annotdir: str, male_threshold: float = 0.8, female_threshold: float = 0.5) -> hl.MatrixTable:
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
    sex_ht = hl.impute_sex(mtx_unphased.GT, aaf_threshold=0.05, female_threshold=female_threshold, male_threshold=male_threshold)
    #export
    sex_ht.export(os.path.join(annotdir, 'mk43_sex_annotated.sex_check.txt.bgz')) # output
    #annotate input (all chroms) mt with imputed sex and write to file
    sex_colnames = ['f_stat', 'is_female']
    sex_ht = sex_ht.select(*sex_colnames)
    mt = mt.annotate_cols(**sex_ht[mt.col_key])
    sex_mt_file = os.path.join(mtdir, "mk43_mt_sex_annotated.mt") # output
    print("Writing to " + sex_mt_file)
    mt.write(sex_mt_file, overwrite=True)

    return mt


def identify_inconsistencies(mt: hl.MatrixTable, mtdir: str, annotdir: str, resourcedir: str):
    '''
    Find samples where annotated sex conflicts with the sex in our metadata
    Find samples where sex is not annotated
    Find samples where f_stat is between 0.2 and 0.8
    :param MatrixTable mt: MT containing imputed sex in column 'sex_check'
    :param str mtdir: directory output matrix tables are written to
    :param str annotdir: directory annotation files are written to
    :param str resourcedir: directory annotation files are written to
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

    metadata_file =  os.path.join(resourcedir, 'mlwh_sample_and_sex.txt')  # resource
    metadata_ht = hl.import_table(metadata_file, delimiter="\t").key_by('accession_number')
    #we only want those from the metadata file where sex is known
    metadata_ht = metadata_ht.filter((metadata_ht.gender == 'Male') | (metadata_ht.gender == 'Female'))

    #annotate the sex-predictions with the manifest sex annotation - need to use a join here
    ht_joined = qc_ht.annotate(manifest_sex = metadata_ht[qc_ht.s].gender)

    #identify samples where imputed sex and manifest sex conflict
    conflicting_sex_ht = ht_joined.filter(((ht_joined.sex == 'male') & (ht_joined.manifest_sex == 'Female')) | (
        (ht_joined.sex == 'female') & (ht_joined.manifest_sex == 'Male')))
    conflicting_sex_ht.export(os.path.join(annotdir, 'mk43_conflicting_sex.txt.bgz')) # output

    #identify samples where f stat is between 0.2 and 0.8
    f_stat_ht = qc_ht.filter( (qc_ht.f_stat > 0.2) & (qc_ht.f_stat < 0.8) )
    f_stat_ht.export(os.path.join(annotdir, 'mk43_sex_annotation_f_stat_outliers.txt.bgz')) # output
    

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

    mt_in_file = os.path.join(mtdir, "mk43_gatk_unprocessed.mt") # input from 1.1
    print("Reading input matrix")
    mt_unfiiltered = hl.read_matrix_table(mt_in_file)

    #apply hard fitlers
    mt_filtered = apply_hard_filters(mt_unfiiltered, mtdir)

    #impute sex
    # TODO: move male_threshold to config?
    mt_sex = impute_sex(mt_filtered, mtdir, annotdir, male_threshold=0.6)

    # annotate_ambiguous_sex(mt_sex, mtdir)
    # TODO: make this optional and check how it affects the downstream steps
    # there is no metadata for our contrived test datasets

    if os.path.exists(resourcedir):
        identify_inconsistencies(mt_sex, mtdir, annotdir, resourcedir)


if __name__ == '__main__':
    main() 

