# apply gnomad's hard filters and impute sex
# input gatk_unprocessed.mt from step 1.1
import argparse
import os

import hail as hl
import pyspark
from utils.utils import parse_config
from wes_qc import hail_utils


def apply_hard_filters(mt: hl.MatrixTable, mtdir: str) -> None:
    """
    Applies hard filters and annotates samples in the filtered set with call rate
    :param MatrixTable mt: MT containing samples to be ascertained for sex
    :param str mtdir: directory output matrix tables are written to
    :return: MatrixTable with hard filtering annotation
    :rtype: MatrixTable
    """
    print("Applying hard filters")
    filtered_mt_file = "file://" + os.path.join(mtdir, "mt_hard_filters_annotated.mt")
    # Keeping only variations that:
    mt = mt.filter_rows(
        (hl.len(mt.alleles) == 2)
        & hl.is_snp(mt.alleles[0], mt.alleles[1])  # Are SNPs
        & (hl.agg.mean(mt.GT.n_alt_alleles()) / 2 > 0.001)  # Have alternate allele in at least 0.1% of samples
        & (hl.agg.fraction(hl.is_defined(mt.GT)) > 0.99)  # Have defined genotype at least in 99% of samples
    )
    mt = mt.annotate_cols(callrate=hl.agg.fraction(hl.is_defined(mt.GT)))
    mt.write(filtered_mt_file, overwrite=True)


def impute_sex(
    mt: hl.MatrixTable, mtdir: str, annotdir: str, male_threshold: float = 0.8, female_threshold: float = 0.5
) -> None:
    """
    Imputes sex, exports data, and annotates mt with this data
    :param MatrixTable mt: MT containing samples to be ascertained for sex
    :param str mtdir: directory output matrix tables are written to
    :param str annotdir: directory annotation files are written to
    :return: MatrixTable with imputed sex annotations stashed in column annotation 'sex_check'
    :rtype: MatrixTable
    """
    print(
        "Imputing sex with male_threshold = " + str(male_threshold) + " and female threshold = " + str(female_threshold)
    )

    # filter to X and select unphased diploid genotypes
    mt1 = hl.filter_intervals(mt, [hl.parse_locus_interval("chrX")])
    mtx_unphased = mt1.select_entries(GT=hl.unphased_diploid_gt_index_call(mt1.GT.n_alt_alleles()))
    # impute sex on the unphased diploid GTs
    sex_ht = hl.impute_sex(
        mtx_unphased.GT, aaf_threshold=0.05, female_threshold=female_threshold, male_threshold=male_threshold
    )
    # export
    sex_ht.export("file://" + os.path.join(annotdir, "sex_annotated.sex_check.csv.gz"))
    # annotate input (all chroms) mt with imputed sex and write to file
    sex_colnames = ["f_stat", "is_female"]
    sex_ht = sex_ht.select(*sex_colnames)
    mt = mt.annotate_cols(**sex_ht[mt.col_key])
    sex_mt_file = "file://" + os.path.join(mtdir, "mt_sex_annotated.mt")
    print("Writing to " + sex_mt_file)
    mt.write(sex_mt_file, overwrite=True)


def identify_inconsistencies(
    mt: hl.MatrixTable, self_reported_sex: str, conflicting_sex: str, f_stat_outliers: str
) -> None:
    """
    Find samples where annotated sex conflicts with the sex in our metadata
    Find samples where sex is not annotated
    Find samples where f_stat is between 0.2 and 0.8
    :param MatrixTable mt: MT containing imputed sex in column 'sex_check'
    :param str mtdir: directory output matrix tables are written to
    :param str annotdir: directory annotation files are written to
    :param str resourcedir: directory annotation files are written to
    """
    print("Annotating samples with inconsistencies:")
    qc_ht = mt.cols()
    # convert is_female boolean to sex
    sex_expr = hl.case().when(qc_ht.is_female, "female").when(qc_ht.is_female == False, "male").default("undetermined")
    qc_ht = qc_ht.annotate(sex=sex_expr).key_by("s")

    # annotate with manifest sex - keyed on ega to match identifiers in matrixtable
    metadata_ht = hl.import_table("file://" + self_reported_sex, delimiter=",").key_by("SampleID")
    # we only want those from the metadata file where sex is known
    metadata_ht = metadata_ht.filter((metadata_ht.Sex == "Male") | (metadata_ht.Sex == "Female"))

    # annotate the sex-predictions with the manifest sex annotation - need to use a join here
    ht_joined = qc_ht.annotate(manifest_sex=metadata_ht[qc_ht.s].Sex)

    # identify samples where imputed sex and manifest sex conflict
    conflicting_sex_ht = ht_joined.filter(
        ((ht_joined.sex == "male") & (ht_joined.manifest_sex == "Female"))
        | ((ht_joined.sex == "female") & (ht_joined.manifest_sex == "Male"))
    )
    conflicting_sex_ht.export("file://" + conflicting_sex)

    # identify samples where f stat is between 0.2 and 0.8
    f_stat_ht = qc_ht.filter((qc_ht.f_stat > 0.2) & (qc_ht.f_stat < 0.8))
    f_stat_ht.export("file://" + f_stat_outliers)


def main() -> None:
    # set up
    inputs = parse_config()

    data_root: str = inputs["data_root"]
    dataset_name: str = inputs["dataset_name"]

    mtdir: str = os.path.join(data_root, inputs["matrixtables_lustre_dir"])
    annot_dir: str = os.path.join(data_root, inputs["annotation_lustre_dir"])

    # initialise hail
    tmp_dir = inputs["tmp_dir"]
    sc = hail_utils.init_hl(tmp_dir)

    # === Determinind state to work ====

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--hard-filter", help="Apply hard filters", action="store_true")
    parser.add_argument("-i", "--impute-sex", help="Impute sex information", action="store_true")
    parser.add_argument("-c", "--check-consistency", help="Annotate and filter merged mt", action="store_true")
    parser.add_argument("-r", "--run", help="Run all steps", action="store_true")
    args = parser.parse_args()

    if args.hard_filter or args.run:
        mt_in_file = "file://" + os.path.join(mtdir, "gatk_unprocessed.mt")
        mt_unfiltered = hl.read_matrix_table(mt_in_file)
        # apply hard fitlers
        apply_hard_filters(mt_unfiltered, mtdir)

    # impute sex
    if args.impute_sex or args.run:
        filtered_mt_file = "file://" + os.path.join(mtdir, "mt_hard_filters_annotated.mt")
        mt_filtered = hl.read_matrix_table(filtered_mt_file)
        impute_sex(mt_filtered, mtdir, annot_dir, male_threshold=0.6)

    # Check for conflicting sex
    if args.check_consistency or args.run:
        self_reported_sex = os.path.join(annot_dir, inputs["self_reported_sex"])
        conflicting_sex = os.path.join(annot_dir, "conflicting_sex.tsv.gz")
        f_stat_outliers = os.path.join(annot_dir, "sex_f_stat_outliers.tsv.gz")

        sex_mt_file = "file://" + os.path.join(mtdir, "mt_sex_annotated.mt")
        mt_sex = hl.read_matrix_table(sex_mt_file)

        identify_inconsistencies(mt_sex, self_reported_sex, conflicting_sex, f_stat_outliers)

    hail_utils.stop_hl(sc)


if __name__ == "__main__":
    main()
