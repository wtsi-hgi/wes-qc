# perform hail sample QC stratified by superpopulation and identify outliers
import hail as hl
from wes_qc import hail_utils
from utils.utils import parse_config, path_spark
from gnomad.sample_qc.filtering import compute_stratified_metrics_filter


# TODO: rename to annotate_with_pop
def annotate_mt(raw_mt_file: str, pop_ht_file: str, annotated_mt_file: str, config: dict = None):
    """
    Annotate mt with superpopulation and sequencing runid
    :param str raw_mt_file: raw mt file
    :param str pop_ht_file: file with population annotations
    :param str runid_file: metadata file to annotate run ids TODO
    :param str annotated_mt_file: annotated mt file
    :param str pop_pandas_file: tsv from pandas df - EGA and pop
    :param dict config: A config object. No effect.

    ### Config fields
    None

    ### Indirect config fields
    step1.gatk_mt_outfile : input path : used in main
    step2.predict_populations.pop_ht_file : input path : used in main
    step2.annotate_with_pop.annotated_mt_file : output path : used in main
    """
    mt = hl.read_matrix_table(path_spark(raw_mt_file))
    pop_ht = hl.read_table(path_spark(pop_ht_file))
    mt = mt.annotate_cols(assigned_pop=pop_ht[mt.s].pop)
    mt.write(path_spark(annotated_mt_file), overwrite=True)  # output
    return mt


def stratified_sample_qc(
    annotated_mt_file: str, mt_qc_outfile: str, ht_qc_cols_outfile: str, qc_filter_file: str, config: dict
):
    """
    Run sample QC and stratify by population
    param str annotated_mt_file: population and run id annotated MT file
    param str mt_qc_outfile: sample QC MT file
    param str ht_qccols_outfile: sample QC columns HT file
    param str qc_filter_file: output file for stratified sample QC HT
    param dict config: config object

    TODO: note about `if min_dp > 0 or min_gq > 0 or min_vaf > 0`

    ### Config fields
    step2.stratified_sample_qc.min_dp : float : TODO
    step2.stratified_sample_qc.min_gq : float : TODO
    step2.stratified_sample_qc.min_vaf : float : TODO
    step2.stratified_sample_qc.output_text_file : output path : TODO
    step2.stratified_sample_qc.output_globals_json_file : output path : TODO

    ### Indirect config fields
    step2.annotate_with_pop.annotated_mt_file : input path : used in main
    step2.stratified_sample_qc.mt_qc_outfile : output path : used in main
    step2.stratified_sample_qc.ht_qc_cols_outfile : output path : used in main
    step2.stratified_sample_qc.qc_filter_file : output path : used in main
    """
    conf = config["step2"]["stratified_sample_qc"]

    mt = hl.read_matrix_table(path_spark(annotated_mt_file))
    # filter to autosomes only
    mt = mt.filter_rows(mt.locus.in_autosome())

    # filter MT by depth/gq/vaf

    min_dp = conf["min_depth"]
    min_gq = conf["min_genotype_quality"]
    min_vaf = conf["min_vaf"]
    if min_dp > 0 or min_gq > 0 or min_vaf > 0:
        vaf = mt.AD[1] / hl.sum(mt.AD)
        print(f"Filtering input MT by depth: DP={min_dp}, genotype quality: GQ={min_gq}, VAF: VAF={min_vaf}")
        filter_condition = (
            (mt.GT.is_het() & (vaf > min_vaf) & (mt.DP > min_dp) & (mt.GQ > min_gq))
            | (mt.GT.is_hom_ref() & (mt.DP > min_dp) & (mt.GQ > min_gq))
            | (mt.GT.is_hom_var() & (mt.DP > min_dp) & (mt.GQ > min_gq))
        )
        fraction_filtered = mt.aggregate_entries(hl.agg.fraction(~filter_condition))
        print(f"Filtering {fraction_filtered * 100:.2f}% entries out of downstream analysis.")
        mt = mt.filter_entries(filter_condition)

    # run sample QC
    mt_with_sampleqc = hl.sample_qc(mt, name="sample_qc")
    # annotate with heterozygosity rate
    mt_with_sampleqc = mt_with_sampleqc.annotate_cols(
        sample_qc=mt_with_sampleqc.sample_qc.annotate(
            heterozygosity_rate=mt_with_sampleqc.sample_qc.n_het / mt_with_sampleqc.sample_qc.n_called
        )
    )

    mt_with_sampleqc.write(path_spark(mt_qc_outfile), overwrite=True)  # output
    mt_with_sampleqc.cols().write(path_spark(ht_qc_cols_outfile), overwrite=True)  # output

    # stratify by pop
    pop_ht = hl.read_table(path_spark(ht_qc_cols_outfile))
    # TODO: shall we move qc_metrics to the config?
    qc_metrics = [
        "heterozygosity_rate",
        "n_snp",
        "r_ti_tv",
        "n_transition",
        "n_transversion",
        "r_insertion_deletion",
        "n_insertion",
        "n_deletion",
        "r_het_hom_var",
    ]

    pop_filter_ht = compute_stratified_metrics_filter(
        pop_ht,
        qc_metrics={metric: pop_ht.sample_qc[metric] for metric in qc_metrics},
        strata={"qc_pop": pop_ht.assigned_pop},
    )

    globals = hl.eval(pop_filter_ht.globals.qc_metrics_stats)

    pop_ht = pop_ht.annotate_globals(qc_metrics_stats=globals)
    pop_ht = pop_ht.annotate(**pop_filter_ht[pop_ht.key]).persist()
    checkpoint = pop_ht.aggregate(hl.agg.count_where(hl.len(pop_ht.qc_metrics_filters) == 0))
    print(f"{checkpoint} exome samples found passing pop filtering")
    pop_ht.write(path_spark(qc_filter_file), overwrite=True)  # output

    # output_text_file = os.path.join(annotdir, "sample_qc_by_pop.tsv.bgz")
    output_text_file = conf["output_text_file"]
    pop_ht.export(path_spark(output_text_file), delimiter="\t")  # output

    # output_globals_json = os.path.join(annotdir, "sample_qc_by_pop.globals.json")
    output_globals_json = conf["output_globals_json_file"]
    pop_ht.globals.export(path_spark(output_globals_json))  # output


def main():
    # set up
    config = parse_config()

    # initialise hail
    tmp_dir = config["general"]["tmp_dir"]
    hail_utils.init_hl(tmp_dir)

    # annotate mt with runid and pop
    raw_mt_file = config["step1"]["gatk_mt_outfile"]
    pop_ht_file = config["step2"]["predict_pops"]["pop_ht_outfile"]
    annotated_mt_file = config["step2"]["annotate_with_pop"]["annotated_mt_file"]
    annotate_mt(raw_mt_file, pop_ht_file, annotated_mt_file)

    # run sample QC and stratify by population
    mt_qc_outfile = config["step2"]["stratified_sample_qc"]["mt_qc_outfile"]
    ht_qc_cols_outfile = config["step2"]["stratified_sample_qc"]["ht_qc_cols_outfile"]
    qc_filter_file = config["step2"]["stratified_sample_qc"]["qc_filter_file"]
    stratified_sample_qc(annotated_mt_file, mt_qc_outfile, ht_qc_cols_outfile, qc_filter_file, config)


if __name__ == "__main__":
    main()
