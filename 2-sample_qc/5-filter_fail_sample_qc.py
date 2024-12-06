# Filter out samples which have failed sample QC
# Filter to remove samples sequenced at Broad
import hail as hl
from utils.utils import parse_config, path_spark
from wes_qc import hail_utils


# TODO: not used? make it optional and mention in cli help
def filter_to_sanger_only(annotated_mt_file: str, sanger_mt_file: str):
    """
    param str annotated_mt_file: File containing MatrixTable annotated with sequencing location
    param str sanger_mt_file: File for sanger only MatrixTable
    population and run id but not filtered
    :return: MatrixTable containing only Sanger samples
    :rtype: hl.MatrixTable
    """
    mt = hl.read_matrix_table(annotated_mt_file)
    mt = mt.filter_cols(mt.sequencing_location == "Sanger")  # filter to Sanger only
    mt.write(sanger_mt_file)


def remove_sample_qc_fails(
    mt_file: str,
    qc_filter_ht_file: str,
    samples_failing_qc_file: str,
    sample_qc_filtered_mt_file: str,
    config: dict = None,
):
    """
    param str mt_file: Input MatrixTable file
    param str qc_filter_ht_file: File contaning sample QC output
    param str samples_failing_qc_file: Output file for list of samples failing QC
    param str sample_qc_filtered_mt_file: Output file for MatrixTable with sample QC fails removed
    param dict config: A config object. No effect.

    ### Config fields
    None

    ### Indirect config fields
    step2.annotate_with_pop.annotated_mt_file : input path : used in main
    step2.stratified_sample_qc.mt_qc_outfile : input path : used in main
    step2.remove_sample_qc_fails.samples_failing_qc_file : output path : used in main
    step2.remove_sample_qc_fails.sample_qc_filtered_mt_file : output path : used in main
    """
    mt = hl.read_matrix_table(path_spark(mt_file))
    sample_qc_ht = hl.read_table(path_spark(qc_filter_ht_file))
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
    fails.export(path_spark(samples_failing_qc_file))  # output
    # filter the sangermt to remove samples that have failed
    fails = fails.key_by("s")
    mt = mt.filter_cols(hl.is_defined(fails[mt.s]), keep=False)
    # save the filtered mt to file
    mt.write(path_spark(sample_qc_filtered_mt_file), overwrite=True)  # output


def main():
    # set up
    config = parse_config()
    tmp_dir = config["general"]["tmp_dir"]

    # = STEP LOGIC = #
    hail_utils.init_hl(tmp_dir)

    # annotated_mt_file = os.path.join(mtdir,"gatk_unprocessed_with_pop.mt")  # annotated but unfiltered mt
    # qc_filter_ht_file = os.path.join(mtdir,"mt_pops_QC_filters.ht")
    # samples_failing_qc_file = os.path.join(annotdir, "samples_failing_qc.tsv.bgz")
    # sample_qc_filtered_mt_file = os.path.join(mtdir,"mt_pops_QC_filters_after_sample_qc.mt")
    annotated_mt_file = config["step2"]["annotate_with_pop"]["annotated_mt_file"]  # annotated but unfiltered mt
    qc_filter_ht_file = config["step2"]["stratified_sample_qc"]["qc_filter_file"]
    samples_failing_qc_file = config["step2"]["remove_sample_qc_fails"]["samples_failing_qc_file"]
    sample_qc_filtered_mt_file = config["step2"]["remove_sample_qc_fails"]["sample_qc_filtered_mt_file"]
    remove_sample_qc_fails(annotated_mt_file, qc_filter_ht_file, samples_failing_qc_file, sample_qc_filtered_mt_file)


if __name__ == "__main__":
    main()
