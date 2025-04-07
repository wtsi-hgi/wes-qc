# apply gnomAD's impute sex
from typing import Any, Optional

import hail as hl
import bokeh
from utils.utils import parse_config, path_spark

from wes_qc import hail_utils


def apply_hard_filters(
    mt: hl.MatrixTable, n_alt_alleles_threshold, defined_gt_frac_threshold, **kwargs
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


def impute_sex(mt: hl.MatrixTable, hail_impute_sex_params: dict[str, Any], **kwargs) -> (hl.MatrixTable, hl.Table):
    """
    Imputes sex, exports data, and annotates mt with this data
    """
    print("===Imputing sex ===")
    mt1 = hl.split_multi_hts(mt)
    mtx_unphased = mt1.select_entries(GT=hl.unphased_diploid_gt_index_call(mt1.GT.n_alt_alleles()))

    # Impute sex on the unphased diploid GTs
    sex_ht = hl.impute_sex(mtx_unphased.GT, **hail_impute_sex_params)

    # convert is_female boolean to sex
    sex_expr = hl.if_else(hl.is_defined(sex_ht.is_female), hl.if_else(sex_ht.is_female, "female", "male"), "undefined")
    sex_ht = sex_ht.annotate(imputed_sex=sex_expr)

    # export
    #

    # annotate input (all chroms) mt with imputed sex and write to file
    sex_colnames = ["f_stat", "is_female", "imputed_sex"]
    # sex_ht = sex_ht.select(*sex_colnames)
    mt = mt.annotate_cols(**sex_ht.select(*sex_colnames)[mt.col_key])

    return mt, sex_ht


def identify_inconsistencies(
    mt: hl.MatrixTable,
) -> Optional[hl.Table]:
    """
    Find samples where annotated sex conflicts with the sex in our metadata
    """
    qc_ht = mt.cols()
    if "self_reported_sex" not in qc_ht.row:
        print("=== No self-reported sex assigned - skipping sex consistency checking ===")
        return None

    print("Annotating samples with inconsistencies:")
    conflicting_sex_ht = qc_ht.filter(
        ((qc_ht.imputed_sex == "male") & (qc_ht.self_reported_sex == "female"))
        | ((qc_ht.imputed_sex == "female") & (qc_ht.self_reported_sex == "male"))
    )
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


def plot_f_stat_histogram(sex_ht: hl.Table, fstat_low: float, fstat_high: float, **kwargs):
    """
    Plot a histogram of F-statistic values from a Hail Table and annotate it with thresholds
    and the count of outliers.

    Parameters:
        sex_ht (hl.Table): Input Hail Table containing the F-statistic data.
        fstat_low (float): Lower F-stat threshold for detecting outliers.
        fstat_high (float): Upper F-stat threshold for detecting outliers.

    Returns:
        bokeh.plotting.figure.Figure: Bokeh figure containing the histogram with thresholds
        and annotation.
    """
    n_outliers = sex_ht.filter((sex_ht.f_stat > fstat_low) & (sex_ht.f_stat < fstat_high)).count()
    plot = hl.plot.histogram(sex_ht.f_stat, legend="Fstat")
    hline_lower = bokeh.models.Span(location=fstat_low, dimension="height", line_color="red", line_width=2)
    hline_higher = bokeh.models.Span(location=fstat_high, dimension="height", line_color="red", line_width=2)
    plot.add_layout(hline_lower)
    plot.add_layout(hline_higher)
    label = bokeh.models.Label(
        x=(fstat_low + fstat_high) / 2, y=10, text=str(n_outliers), text_font_size="12pt", text_color="red"
    )  # Add label to the figure
    plot.add_layout(label)
    return plot


def main():
    ## = STEP SETUP = #
    config = parse_config()
    tmp_dir = config["general"]["tmp_dir"]

    # = STEP DEPENDENCIES = #
    mt_infile = config["step1"]["validate_gtcheck"]["mt_gtcheck_validated"]  # input from 1.3

    # = STEP OUTPUTS = #
    sex_mt_file = config["step2"]["impute_sex"]["sex_mt_outfile"]
    sex_ht_outfile = config["step2"]["impute_sex"]["sex_ht_outfile"]
    conflicting_sex_report_file = config["step2"]["sex_inconsistencies"]["conflicting_sex_report_file"]
    fstat_outliers_report_file = config["step2"]["f_stat_outliers"]["fstat_outliers_report_file"]
    fstat_hist_path = config["step2"]["f_stat_outliers"]["fstat_hist_path"]

    # = STEP LOGIC = #
    _ = hail_utils.init_hl(tmp_dir)

    print("Reading input matrix")
    mt_unfiltered = hl.read_matrix_table(path_spark(mt_infile))

    # apply hard filters
    mt_filtered = apply_hard_filters(mt_unfiltered, **config["step2"]["sex_annotation_hard_filters"])

    # impute sex
    mt_sex, sex_ht = impute_sex(mt_filtered, **config["step2"]["impute_sex"])
    sex_ht.export(path_spark(sex_ht_outfile))

    # Save f-stat outliers in the separate file
    f_stat_ht_outliers = select_fstat_outliers(sex_ht, **config["step2"]["f_stat_outliers"])
    f_stat_ht_outliers.export(path_spark(fstat_outliers_report_file))  # output

    # Plot F-stat histogram
    fstat_hist = plot_f_stat_histogram(sex_ht, **config["step2"]["f_stat_outliers"])
    bokeh.io.output_file(fstat_hist_path)
    bokeh.io.save(fstat_hist)

    print("--- Writing to " + sex_mt_file)
    mt_sex.write(path_spark(sex_mt_file), overwrite=True)

    # identify and report inconsistencies
    conflicting_sex_ht = identify_inconsistencies(mt_sex)
    if conflicting_sex_ht is not None:
        conflicting_sex_ht.export(path_spark(conflicting_sex_report_file))  # output


if __name__ == "__main__":
    main()
