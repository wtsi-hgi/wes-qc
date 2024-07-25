# Filter out samples which have failed sample QC
# Filter to remove samples sequenced at Broad
import os

import hail as hl
import pyspark
from utils.utils import parse_config
from wes_qc import hail_utils


def filter_to_sanger_only(annotated_mt_file: str, sanger_mt_file: str) -> None:
    """
    param str annotated_mt_file: File containing MatrixTable annotated with sequencing location
    param str sanger_mt_file: File for sanger only MatrixTable
    population and run id but not filtered
    :return: MatrixTable containing only Sanger samples
    :rtype: hl.MatrixTable
    """
    mt = hl.read_matrix_table(annotated_mt_file)
    mt = mt.filter_cols(mt.sequencing_location == "Sanger")  # filter to Sanger only
    mt.write(sanger_mt_file, overwrite=True)


def remove_sample_qc_fails(
    annotated_mt_file: str,
    qc_filter_ht_file: str,
    verifybamid_file: str,
    samples_failing_qc_file: str,
    sample_qc_filtered_mt_file: str,
    missed_freemix_file: str,
) -> None:
    """
    param str sanger_mt_file: Input MatrixTable file
    param str qc_filter_ht_file: File contaning sample QC output
    param str samples_failing_qc_file: Output file for list of samples failing QC
    param str verifybamid_file: VerfiyBamID2 results for all samples
    param str gtcheck_fail_file: Samples which fail on gtcheck
    param str sample_qc_filtered_mt_file: Output file for MatrixTable with sample QC fails removed
    """
    annotated_mt = hl.read_matrix_table(annotated_mt_file)
    n_samples_initial = annotated_mt.count_cols()
    print(f"=== Loaded: {n_samples_initial} samples")
    sample_qc_ht = hl.read_table(qc_filter_ht_file)
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
    fails.export(samples_failing_qc_file)
    # filter the annotated_mt to remove samples that have failed
    fails = fails.key_by("s")
    annotated_mt = annotated_mt.filter_cols(hl.is_defined(fails[annotated_mt.s]), keep=False)
    print(f"=== Survived after QC filtering: {annotated_mt.count_cols()} samples")

    # remove samples failing FREEMIX check
    verifybamid_ht = hl.import_table(verifybamid_file, impute=True)
    verifybamid_ht = verifybamid_ht.transmute(s=verifybamid_ht["#SEQ_ID"])
    verifybamid_ht = verifybamid_ht.key_by("s")
    annotated_mt = annotated_mt.annotate_cols(freemix=verifybamid_ht[annotated_mt.s]["FREEMIX"])

    sample_table = annotated_mt.cols()
    missed_freemix = sample_table.filter(hl.is_defined(sample_table.freemix), keep=False)
    if missed_freemix.count() > 0:
        print(f"=== WARNING: Found samples without freemix scores: {missed_freemix_file}")
        missed_freemix.export(missed_freemix_file)
    else:
        print("=== All samples have freemix scores")
    annotated_mt = annotated_mt.filter_cols(hl.is_defined(annotated_mt.freemix))
    annotated_mt = annotated_mt.filter_cols(annotated_mt.freemix <= 0.05)

    # save the filtered mt to file
    annotated_mt.write(sample_qc_filtered_mt_file, overwrite=True)
    n_samples_final = annotated_mt.count_cols()
    print(f"=== Survived after all filtering: {n_samples_final} samples")


def main() -> None:
    # set up
    inputs = parse_config()

    data_root: str = inputs["data_root"]
    dataset_name: str = inputs["dataset_name"]
    mtdir: str = os.path.join(data_root, inputs["matrixtables_lustre_dir"])
    resourcesdir = os.path.join(data_root, inputs["resource_dir"])
    annotdir: str = os.path.join(data_root, inputs["annotation_lustre_dir"])
    plotdir = os.path.join(data_root, inputs["plots_lustre_dir"])

    # initialise hail
    sc = hail_utils.init_hl(inputs["tmp_dir"])

    annotated_mt_file = os.path.join(mtdir, "gatk_unprocessed_with_pop.mt")  # annotated but unfiltered mt
    qc_filter_ht_file = os.path.join(mtdir, "mt_pops_QC_filters.ht")
    verifybamid_file = os.path.join(annotdir, "verify_bam_id_result_concat.selfSM")

    sample_qc_filtered_mt_file = os.path.join(mtdir, "mt_pops_QC_filters_after_sample_qc.mt")

    samples_failing_qc_file = os.path.join(annotdir, "samples_failing_qc.tsv.gz")
    missed_freemix_file = os.path.join(annotdir, "samples_missed_freemix.tsv")

    remove_sample_qc_fails(
        "file://" + annotated_mt_file,
        "file://" + qc_filter_ht_file,
        "file://" + verifybamid_file,
        "file://" + samples_failing_qc_file,
        "file://" + sample_qc_filtered_mt_file,
        "file://" + missed_freemix_file,
    )

    # The previous version with removing samples after failed gtcheck
    # remove_sample_qc_fails(sanger_mt_file, qc_filter_ht_file, samples_failing_qc_file, verifybamid_file, gtcheck_fail_file, sample_qc_filtered_mt_file)

    hail_utils.stop_hl(sc)


if __name__ == "__main__":
    main()
