# Filter out samples which have failed sample QC
# Filter to remove samples sequenced at Broad
import hail as hl
from utils.utils import parse_config, path_spark
from wes_qc import hail_utils


def remove_sample_qc_fails(
    mt: hl.MatrixTable,
    sample_qc_ht: hl.Table,
) -> tuple[hl.MatrixTable, hl.Table]:
    """
    Removes samples from the MatrixTable that fail quality control (QC) checks based on a sample QC Hail Table.
    Additionally, returns a Table of failed samples.

    Args:
        mt (hl.MatrixTable): Input MatrixTable containing genomic data.
        sample_qc_ht (hl.Table): Table containing sample QC metrics and filters used for
            determining sample failures.

    Returns:
        hl.MatrixTable: MatrixTable with samples failing QC removed.
        hl.Table: Table of samples that failed QC.
    """

    # identify samples which have failed any of the metrics tested
    sample_qc_ht = sample_qc_ht.annotate(filter_fail_count=(hl.len(sample_qc_ht.qc_metrics_filters)))

    filter_expr = (
        hl.case()
        .when(sample_qc_ht.filter_fail_count == 0, "Pass")
        .when(sample_qc_ht.filter_fail_count > 0, "Fail")
        .default("")
    )
    sample_qc_ht = sample_qc_ht.annotate(filter_result=filter_expr).key_by("s")
    # save a list of samples that have failed QC
    fail_ht = sample_qc_ht.filter(sample_qc_ht.filter_result == "Fail").key_by()
    fails = fail_ht.select(fail_ht.s)
    fails = fails.key_by("s")

    # Removing samples that are present in the fails table
    mt = mt.filter_cols(hl.is_defined(fails[mt.s]), keep=False)
    return mt, fails


def remove_samples_explicit(mt: hl.MatrixTable, samples_to_remove: hl.Table) -> hl.MatrixTable:
    """
    Filters out specific samples from a given Hail MatrixTable based on a provided Table.

    Args:
        mt (hl.MatrixTable): The MatrixTable from which samples need to be removed.
        samples_to_remove (hl.Table): A Table containing the samples to be removed
            from the MatrixTable. This Table should have keys matching with the
            column key of the MatrixTable.

    Returns:
        hl.MatrixTable: A new MatrixTable with the specified samples removed.
    """
    mt = mt.filter_cols(hl.is_missing(samples_to_remove[mt.s]))
    return mt


def main():
    ## = STEP SETUP = #
    config = parse_config()
    tmp_dir = config["general"]["tmp_dir"]

    # = STEP DEPENDENCIES = #
    annotated_mt_file = config["step2"]["annotate_with_pop"]["annotated_mt_file"]  # annotated but unfiltered mt
    qc_filter_ht_file = config["step2"]["stratified_sample_qc"]["qc_filter_file"]  # Results of stratified filtering
    samples_to_remove_file = config["step2"]["sample_qc_combine_results"]["samples_to_remove_file"]

    # = STEP OUTPUTS = #
    samples_failing_qc_file = config["step2"]["remove_sample_qc_fails"][
        "samples_failing_qc_file"
    ]  # Table of samples failed Stratified QC
    sample_qc_filtered_mt_file = config["step2"]["remove_sample_qc_fails"]["sample_qc_filtered_mt_file"]

    # = STEP LOGIC = #
    hail_utils.init_hl(tmp_dir)

    mt = hl.read_matrix_table(path_spark(annotated_mt_file))
    sample_qc_ht = hl.read_table(path_spark(qc_filter_ht_file))
    print("=== Removing samples failed stratified sample QC ===")
    mt, fails = remove_sample_qc_fails(mt, sample_qc_ht)
    fails.export(path_spark(samples_failing_qc_file))  # output

    if samples_to_remove_file is not None:
        samples_to_remove = (
            hl.import_table(path_spark(samples_to_remove_file), no_header=True, impute=False)
            .rename({"f0": "s"})
            .key_by("s")
        )
        print(f"=== Removing {samples_to_remove.count()} samples listed in {samples_to_remove_file} ===")
        mt = remove_samples_explicit(mt, samples_to_remove)

    mt.write(path_spark(sample_qc_filtered_mt_file), overwrite=True)  # output


if __name__ == "__main__":
    main()
