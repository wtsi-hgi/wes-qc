# apply gnomad's hard filters and impute sex
# input gatk_unprocessed.mt from step 1.1
import hail as hl
from utils.utils import parse_config, path_spark

from wes_qc import hail_utils


def apply_hard_filters(mt: hl.MatrixTable, config: dict) -> hl.MatrixTable:
    """
    Applies hard filters and annotates samples in the filtered set with call rate
    :param MatrixTable mt: MT containing samples to be ascertained for sex
    :param str mtdir: directory output matrix tables are written to
    :param dict config:
    :return: MatrixTable with hard filtering annotation
    :rtype: MatrixTable

    ### Config fields
    step2.sex_annotation_hard_filters.filtered_mt_outfile : path
    step2.sex_annotation_hard_filters.n_alt_alleles_threshold : float
    step2.sex_annotation_hard_filters.defined_gt_frac_threshold : float
    """
    conf = config["step2"]["sex_annotation_hard_filters"]

    print("Applying hard filters")
    filtered_mt_file = path_spark(conf["filtered_mt_outfile"])  # output

    mt = mt.filter_rows(
        (hl.len(mt.alleles) == 2)
        & hl.is_snp(mt.alleles[0], mt.alleles[1])
        & (hl.agg.mean(mt.GT.n_alt_alleles()) / 2 > conf["n_alt_alleles_threshold"])
        & (hl.agg.fraction(hl.is_defined(mt.GT)) > conf["defined_gt_frac_threshold"])
    )
    mt = mt.annotate_cols(callrate=hl.agg.fraction(hl.is_defined(mt.GT)))
    mt.write(filtered_mt_file, overwrite=True)

    return mt


def impute_sex(
    mt: hl.MatrixTable, aaf_threshold: float, female_threshold: float, male_threshold: float, sex_ht_outfile, **kwargs
) -> hl.MatrixTable:
    """
    Imputes sex, exports data, and annotates mt with this data
    :param MatrixTable mt: MT containing samples to be ascertained for sex
    :param dict config:
    :return: MatrixTable with imputed sex annotations stashed in column annotation 'sex_check'
    :rtype: MatrixTable

    ### Config fields
    step2.impute_sex.sex_ht_outfile : path
    step2.impute_sex.sex_mt_outfile : path
    step2.impute_sex.female_threshold : float
    step2.impute_sex.male_threshold : float
    step2.impute_sex.aaf_threshold : float
    """

    print(f"Imputing sex with male_threshold = {male_threshold} and female threshold = {female_threshold}")

    # filter to X and select unphased diploid genotypes - no need to filter to X as impute_sex takes care of this
    # mt1 = hl.filter_intervals(mt, [hl.parse_locus_interval('chrX')])
    mt1 = hl.split_multi_hts(mt)
    mtx_unphased = mt1.select_entries(GT=hl.unphased_diploid_gt_index_call(mt1.GT.n_alt_alleles()))
    # imput sex on the unphased diploid GTs
    sex_ht = hl.impute_sex(
        mtx_unphased.GT,
        aaf_threshold=aaf_threshold,
        female_threshold=female_threshold,
        male_threshold=male_threshold,
    )

    # convert is_female boolean to sex
    sex_expr = hl.if_else(hl.is_defined(sex_ht.is_female), hl.if_else(sex_ht.is_female, "female", "male"), "undefined")
    # sex_expr = hl.case().when(mt.is_female, "female").when(~mt.is_female, "male").default("undetermined")
    sex_ht = sex_ht.annotate(imputed_sex=sex_expr)

    # export
    sex_ht.export(path_spark(sex_ht_outfile))  # output

    # annotate input (all chroms) mt with imputed sex and write to file
    sex_colnames = ["f_stat", "is_female", "imputed_sex"]
    sex_ht = sex_ht.select(*sex_colnames)
    mt = mt.annotate_cols(**sex_ht[mt.col_key])

    # sex_mt_file = path_spark(conf["sex_mt_outfile"])  # output
    # print("Writing to " + sex_mt_file)
    # mt.write(sex_mt_file, overwrite=True)

    return mt


def identify_inconsistencies(mt: hl.MatrixTable, config: dict):
    """
    Find samples where annotated sex conflicts with the sex in our metadata
    Find samples where sex is not annotated
    Find samples where f_stat is between fstat_low and fstat_high
    :param MatrixTable mt: MT containing imputed sex in column 'sex_check'
    :param dict config:

    ### Config fields
    step2.sex_inconsistencies.sex_metadata_file : input path : TODO explain metadata structure and constants
    step2.sex_inconsistencies.conflicting_sex_report_file : output path : TODO
    step2.sex_inconsistencies.fstat_outliers_report_file : output path : TODO
    step2.sex_inconsistencies.fstat_low : float
    step2.sex_inconsistencies.fstat_high : float
    """
    conf = config["step2"]["sex_inconsistencies"]

    qc_ht = mt.cols()
    if "self_reported_sex" not in qc_ht.row:
        print("=== No self-reported sex assined - skipping sex consistency checking ===")
        return

    print("Annotating samples with inconsistencies:")

    # annotate with manifest sex - keyed on ega to match identifiers in matrixtable

    # annotate the sex-predictions with the manifest sex annotation - need to use a join here
    ht_joined = qc_ht

    # identify samples where imputed sex and manifest sex conflict
    conflicting_sex_ht = ht_joined.filter(
        ((ht_joined.imputed_sex == "male") & (ht_joined.self_reported_sex == "female"))
        | ((ht_joined.imputed_sex == "female") & (ht_joined.self_reported_sex == "male"))
    )

    conflicting_sex_ht.export(path_spark(conf["conflicting_sex_report_file"]))  # output

    # identify samples where f stat is between fstat_low and fstat_high
    f_stat_ht = qc_ht.filter((qc_ht.f_stat > conf["fstat_low"]) & (qc_ht.f_stat < conf["fstat_high"]))
    f_stat_ht.export(path_spark(conf["fstat_outliers_report_file"]))  # output


def main():
    ## = STEP SETUP = #
    config = parse_config()
    tmp_dir = config["general"]["tmp_dir"]

    # = STEP DEPENDENCIES = #
    mt_infile = config["step1"]["validate_gtcheck"]["mt_gtcheck_validated"]  # input from 1.3

    # = STEP LOGIC = #
    _ = hail_utils.init_hl(tmp_dir)

    print("Reading input matrix")
    mt_unfiltered = hl.read_matrix_table(path_spark(mt_infile))

    # apply hard fitlers
    mt_filtered = apply_hard_filters(mt_unfiltered, config)

    # impute sex
    mt_sex = impute_sex(mt_filtered, **config["step2"]["impute_sex"])

    # TODO: make this optional and check how it affects the downstream steps
    # there is no metadata for our contrived test datasets
    # identify_inconsistencies
    identify_inconsistencies(mt_sex, config)


if __name__ == "__main__":
    main()
