#apply gnomad's hard filters and impute sex
#input gatk_unprocessed.mt from step 1.1
import os
import hail as hl
import pyspark
    from utils.utils import parse_config, path_local, path_spark, __expand_cvars
import os

def apply_hard_filters(mt: hl.MatrixTable, config: dict) -> hl.MatrixTable:
    '''
    Applies hard filters and annotates samples in the filtered set with call rate
    :param MatrixTable mt: MT containing samples to be ascertained for sex
    :param str mtdir: directory output matrix tables are written to
    :param dict config:
    :return: MatrixTable with hard filtering annotation
    :rtype: MatrixTable

    ### Config fields
    ```
    step2.sex_annotation_hard_filters.filtered_mt_outfile : path
    step2.sex_annotation_hard_filters.n_alt_alleles_threshold : float
    step2.sex_annotation_hard_filters.defined_gt_frac_threshold : float
    ```
    '''
    conf = config['step2']['sex_annotation_hard_filters']

    print("Applying hard filters")
    filtered_mt_file = path_spark(conf['filtered_mt_outfile']) # output

    mt = mt.filter_rows((hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1]) &
        (hl.agg.mean(mt.GT.n_alt_alleles()) / 2 > conf['n_alt_alleles_threshold']) &
        (hl.agg.fraction(hl.is_defined(mt.GT)) > conf['defined_gt_frac_threshold']))
    mt = mt.annotate_cols(callrate=hl.agg.fraction(hl.is_defined(mt.GT)))
    mt.write(filtered_mt_file, overwrite=True)

    return mt


def impute_sex(mt: hl.MatrixTable, config: dict) -> hl.MatrixTable:
    '''
    Imputes sex, exports data, and annotates mt with this data
    :param MatrixTable mt: MT containing samples to be ascertained for sex
    :param dict config:
    :return: MatrixTable with imputed sex annotations stashed in column annotation 'sex_check'
    :rtype: MatrixTable

    ### Config fields
    ```
    step2.impute_sex.sex_ht_outfile : path
    step2.impute_sex.sex_mt_outfile : path
    step2.impute_sex.female_threshold : float
    step2.impute_sex.male_threshold : float
    step2.impute_sex.aaf_threshold : float
    ```
    '''

    conf = config['step2']['impute_sex']
    print("Imputing sex with male_threshold = " + str(conf['male_threshold']) + " and female threshold = " + str(conf['female_threshold']))

    #filter to X and select unphased diploid genotypes - no need to filter to X as impute_sex takes care of this
    #mt1 = hl.filter_intervals(mt, [hl.parse_locus_interval('chrX')])
    mt1 = hl.split_multi_hts(mt)
    mtx_unphased = mt1.select_entries(GT=hl.unphased_diploid_gt_index_call(mt1.GT.n_alt_alleles()))
    #imput sex on the unphased diploid GTs
    sex_ht = hl.impute_sex(mtx_unphased.GT, aaf_threshold=conf['aaf_threshold'], female_threshold=conf['female_threshold'], male_threshold=conf['male_threshold'])
    #export
    sex_ht.export(path_spark(conf['sex_ht_outfile'])) # output
    #annotate input (all chroms) mt with imputed sex and write to file
    sex_colnames = ['f_stat', 'is_female']
    sex_ht = sex_ht.select(*sex_colnames)
    mt = mt.annotate_cols(**sex_ht[mt.col_key])
    sex_mt_file = path_spark(conf['sex_mt_outfile']) # output
    print("Writing to " + sex_mt_file)
    mt.write(sex_mt_file, overwrite=True)

    return mt


def identify_inconsistencies(mt: hl.MatrixTable, config: dict):
    '''
    Find samples where annotated sex conflicts with the sex in our metadata
    Find samples where sex is not annotated
    Find samples where f_stat is between fstat_low and fstat_high
    :param MatrixTable mt: MT containing imputed sex in column 'sex_check'
    :param dict config:

    ### Config fields
    ```
    step2.sex_inconsistencies.sex_metadata_file : input path : TODO explain metadata structure and constants
    step2.sex_inconsistencies.conflicting_sex_report_path : output path : TODO
    step2.sex_inconsistencies.fstat_outliers_report_path : output path : TODO
    step2.sex_inconsistencies.fstat_low : float
    step2.sex_inconsistencies.fstat_high : float
    ```
    '''
    conf = config['step2']['sex_inconsistencies']

    # TODO: do we need such a detailed logging, or a single if (... and ... and ...) will suffice?
    error = False
    if not os.path.exists(conf['sex_metadata_file']):
        print("error: identify_inconsistencies: sex_metadata_file missing")
        error = True
    if not os.path.exists(conf['conflicting_sex_report_path']):
        print("error: identify_inconsistencies: conflicting_sex_report_path missing")
        error = True
    if not os.path.exists(conf['fstat_outliers_report_path']):
        print("error: identify_inconsistencies: fstat_outliers_report_path missing")
        error = True
    if error:
        print("skip identify_inconsistencies because of previous errors")
        return

    print("Annotating samples with inconsistencies:")
    qc_ht = mt.cols()
    #convert is_female boolean to sex
    sex_expr = (hl.case()
            .when(qc_ht.is_female, "female")
            .when(qc_ht.is_female == False, "male")
            .default("undetermined"))
    qc_ht = qc_ht.annotate(sex=sex_expr).key_by('s')

    #annotate with manifest sex - keyed on ega to match identifiers in matrixtable

    metadata_ht = hl.import_table(path_spark(conf['sex_metadata_file']), delimiter="\t").key_by('accession_number')
    #we only want those from the metadata file where sex is known
    metadata_ht = metadata_ht.filter((metadata_ht.gender == 'Male') | (metadata_ht.gender == 'Female'))

    #annotate the sex-predictions with the manifest sex annotation - need to use a join here
    ht_joined = qc_ht.annotate(manifest_sex = metadata_ht[qc_ht.s].gender)

    #identify samples where imputed sex and manifest sex conflict
    conflicting_sex_ht = ht_joined.filter(((ht_joined.sex == 'male') & (ht_joined.manifest_sex == 'Female')) | (
        (ht_joined.sex == 'female') & (ht_joined.manifest_sex == 'Male')))
    conflicting_sex_ht.export(path_spark(conf['conflicting_sex_report_path'])) # output

    #identify samples where f stat is between fstat_low and fstat_high
    f_stat_ht = qc_ht.filter( (qc_ht.f_stat > conf['fstat_low']) & (qc_ht.f_stat < conf['fstat_high']) )
    f_stat_ht.export(path_spark(conf['fstat_outliers_report_path'])) # output
    

def main():
    #set up
    config = parse_config()
    #importmtdir = inputs['load_matrixtables_lustre_dir']

    #initialise hailS
    tmp_dir = config['general']['tmp_dir']
    # sc = pyspark.SparkContext()
    sc = pyspark.SparkContext.getOrCreate()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38", idempotent=True)

    mt_infile = config['step1']['gatk_mt_outfile'] # input from 1.1
    print("Reading input matrix")
    mt_unfiltered = hl.read_matrix_table(mt_infile)

    #apply hard fitlers
    mt_filtered = apply_hard_filters(mt_unfiltered, config)

    #impute sex
    mt_sex = impute_sex(mt_filtered, config)

    # TODO: where is this function?
    # annotate_ambiguous_sex(mt_sex, mtdir)

    # TODO: make this optional and check how it affects the downstream steps
    # there is no metadata for our contrived test datasets
    #identify_inconsistencies
    identify_inconsistencies(mt_sex, config)
    
if __name__ == '__main__':
    main() 

