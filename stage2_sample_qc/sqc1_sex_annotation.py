# apply gnomAD's impute sex
from typing import Any, Optional

import hail as hl
import bokeh
from wxs_qc.config import get_config
from wxs_qc.hail_utils import path_spark

from wxs_qc import hail_utils, constants


def apply_hard_filters(
    mt: hl.MatrixTable, n_alt_alleles_threshold: float, defined_gt_frac_threshold: float, **kwargs
) -> hl.MatrixTable:
    """
    Applies hard filters and annotates samples in the filtered set with call rate
    """
    print("=== Applying hard filters before sex prediction ===")
    mt = mt.filter_rows(
        (hl.len(mt.alleles) == 2)
        & hl.is_snp(mt.alleles[0], mt.alleles[1])
        & (hl.agg.mean(mt.GT.n_alt_alleles()) / 2 > n_alt_alleles_threshold)
        & (hl.agg.fraction(hl.is_defined(mt.GT)) > defined_gt_frac_threshold)
    )
    mt = mt.annotate_cols(callrate=hl.agg.fraction(hl.is_defined(mt.GT)))
    return mt


def impute_sex(
    mt: hl.MatrixTable, n_partitions: int, hail_impute_sex_params: dict[str, Any], **kwargs
) -> tuple[hl.MatrixTable, hl.Table]:
    """
    Imputes sex, exports data, and annotates mt with this data.

    Output contract:
        Materialized sex Table.
        The returned MatrixTable remains lazy.
    """
    print("--- Filtering matrixtable for sex imputation ---")
    mtx_impute = mt.select_entries(GT=mt.GT)

    mtx_unphased = mtx_impute.select_entries(GT=hl.unphased_diploid_gt_index_call(mtx_impute.GT.n_alt_alleles()))
    # Coalesce to a smaller partition number after a heavy filtering
    mtx_unphased = mtx_unphased.repartition(n_partitions, shuffle=False)
    mtx_unphased = mtx_unphased.checkpoint(hl.utils.new_temp_file("mtx_impute_sex_unphased", "mt"), overwrite=True)

    # Impute sex on the unphased diploid GTs
    print("--- Imputing sex ---")
    sex_ht = hl.impute_sex(mtx_unphased.GT, **hail_impute_sex_params)

    # convert is_female boolean to sex
    sex_expr = hl.if_else(hl.is_defined(sex_ht.is_female), hl.if_else(sex_ht.is_female, "female", "male"), "undefined")
    sex_ht = sex_ht.annotate(imputed_sex=sex_expr)
    sex_ht = sex_ht.checkpoint(hl.utils.new_temp_file("sqc1_sex_annotation", "ht"), overwrite=True)

    # annotate input (all chroms) mt with imputed sex and write to file
    sex_colnames = ["f_stat", "is_female", "imputed_sex"]
    # sex_ht = sex_ht.select(*sex_colnames)
    mt = mt.annotate_cols(**sex_ht.select(*sex_colnames)[mt.col_key])

    return mt, sex_ht


def identify_inconsistencies(
    mt: hl.MatrixTable,
) -> Optional[hl.Table]:
    """
    Find samples where imputed sex does not exactly match the self-reported sex.
    Samples with undefined or missing sex are reported as inconsistent.
    """
    qc_ht = mt.cols()
    if "self_reported_sex" not in qc_ht.row:
        print("=== No self-reported sex assigned - skipping sex consistency checking ===")
        return None

    print("=== Annotating samples with inconsistencies: ===")
    # Only an exact male/female match is consistent.
    # Undefined or missing self-reported or imputed sex is reported as inconsistent.
    is_consistent = hl.coalesce(
        (qc_ht.imputed_sex == qc_ht.self_reported_sex) & (qc_ht.imputed_sex != "undefined"),
        False,
    )
    conflicting_sex_ht = qc_ht.filter(~is_consistent)
    return conflicting_sex_ht


def select_fstat_outliers(sex_ht: hl.Table, fstat_low: float, fstat_high: float, **kwargs) -> hl.Table:
    """
    Identify F-stat outliers.

    Arguments:
        sex_ht (hl.Table): Hail Table containing sample information with an 'f_stat'
            field representing the F-stat values for the samples.
        fstat_low (float): Lower bound for F-stat filtering (exclusive).
        fstat_high (float): Upper bound for F-stat filtering (exclusive).

    Returns:
        hl.Table: A filtered Hail Table containing rows where the 'f_stat' lies between
        the specified lower and upper bounds.
    """
    return sex_ht.filter((sex_ht.f_stat > fstat_low) & (sex_ht.f_stat < fstat_high))


def plot_f_stat_histogram(
    sex_ht: hl.Table, fstat_low: float, fstat_high: float, text_size=constants.plots_text_size, **kwargs
):
    """
    Plot a histogram of F-statistic values from a Hail Table and annotate it with thresholds
    and the count of outliers.
    """
    n_outliers = select_fstat_outliers(sex_ht, fstat_low, fstat_high).count()
    plot = hl.plot.histogram(sex_ht.f_stat, legend="Fstat", title="F-stat distribution")
    hline_lower = bokeh.models.Span(location=fstat_low, dimension="height", line_color="red", line_width=2)
    hline_higher = bokeh.models.Span(location=fstat_high, dimension="height", line_color="red", line_width=2)
    plot.add_layout(hline_lower)
    plot.add_layout(hline_higher)
    label = bokeh.models.Label(
        x=(fstat_low + fstat_high) / 2, y=10, text=str(n_outliers), text_font_size="12pt", text_color="red"
    )  # Add label to the figure
    plot.add_layout(label)
    plot.axis.axis_label_text_font_size = text_size
    plot.axis.major_label_text_font_size = text_size
    plot.title.text_font_size = text_size
    return plot


def main():
    ## = STEP SETUP = #
    config = get_config()
    tmp_dir = config["general"]["tmp_dir"]

    # = STEP PARAMETERS = #
    n_partitions = config["general"]["n_partitions"]

    # = STEP DEPENDENCIES = #
    mt_infile = config["stage1"]["validate_gtcheck"]["mt_gtcheck_validated"]  # input from 1.3

    # = STEP OUTPUTS = #
    sex_mt_file = config["stage2"]["impute_sex"]["sex_mt_outfile"]
    sex_ht_outfile = config["stage2"]["impute_sex"]["sex_ht_outfile"]
    conflicting_sex_report_file = config["stage2"]["sex_inconsistencies"]["conflicting_sex_report_file"]
    fstat_outliers_report_file = config["stage2"]["f_stat_outliers"]["fstat_outliers_report_file"]
    fstat_hist_path = config["stage2"]["f_stat_outliers"]["fstat_hist_path"]

    # = STEP LOGIC = #
    _ = hail_utils.init_hl(tmp_dir)

    print("Reading input matrix")
    mt_unfiltered = hl.read_matrix_table(path_spark(mt_infile))

    # apply hard filters
    mt_filtered = apply_hard_filters(mt_unfiltered, **config["stage2"]["sex_annotation_hard_filters"])

    # impute sex
    mt_sex, sex_ht = impute_sex(mt_filtered, n_partitions, **config["stage2"]["impute_sex"])
    print("--- Writing to " + sex_mt_file)
    mt_sex.write(path_spark(sex_mt_file), overwrite=True)
    sex_ht.export(path_spark(sex_ht_outfile))

    # Save f-stat outliers in the separate file
    f_stat_ht_outliers = select_fstat_outliers(sex_ht, **config["stage2"]["f_stat_outliers"])
    f_stat_ht_outliers.export(path_spark(fstat_outliers_report_file))  # output

    # Plot F-stat histogram
    fstat_hist = plot_f_stat_histogram(sex_ht, **config["stage2"]["f_stat_outliers"])
    bokeh.io.output_file(fstat_hist_path)
    bokeh.io.save(fstat_hist)

    # identify and report inconsistencies
    conflicting_sex_ht = identify_inconsistencies(mt_sex)
    if conflicting_sex_ht is not None:
        conflicting_sex_ht.export(path_spark(conflicting_sex_report_file))  # output


if __name__ == "__main__":
    main()
