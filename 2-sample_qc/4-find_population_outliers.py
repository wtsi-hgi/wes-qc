# perform hail sample QC stratified by superpopulation and identify outliers
import os

import hail as hl
import pyspark
from utils.utils import parse_config
from gnomad.sample_qc.filtering import compute_stratified_metrics_filter

from wes_qc import hail_utils


def annotate_mt(raw_mt_file: str, pop_ht_file: str, annotated_mt_file: str) -> None:
    """
    Annotate mt with superpopulation and sequencing runid
    :param str raw_mt_file: raw mt file
    :param str pop_ht_file: file with population annotations
    :param str runid_file: metadata file to annotate run ids
    :param str annotated_mt_file: annotated mt file
    """
    mt = hl.read_matrix_table(raw_mt_file)
    pop_ht = hl.read_table(pop_ht_file)
    mt = mt.annotate_cols(assigned_pop=pop_ht[mt.s].pop)
    mt.write(annotated_mt_file, overwrite=True)


def stratified_sample_qc(
    annotated_mt_file: str, mt_qc_outfile: str, ht_qc_cols_outfile: str, qc_filter_file: str, sampleqc_report: str
) -> None:
    """
    Run sample QC and stratify by population
    param str annotated_mt_file: population and run id annotated MT file
    param str mt_qc_outfile: sample QC MT file
    param str ht_qccols_outfile: sample QC columns HT file
    param str qc_filter_file: output file for stratified sample QC HT
    param str annotdir: output directory for annotations
    """
    mt = hl.read_matrix_table(annotated_mt_file)
    # filter to autosomes only
    mt = mt.filter_rows(mt.locus.in_autosome())
    n_samples = mt.count_cols()
    print(f"{n_samples} exome samples before filtering")

    # filter MT by depth/gq/vaf
    # We do this filtering so that samples are not penalised by the presence of low quality variants.
    # Usially there is no need to change these filters
    min_dp = 20
    min_gq = 20
    min_vaf = 0.25
    if min_dp > 0 or min_gq > 0 or min_vaf > 0:
        vaf = mt.AD[1] / hl.sum(mt.AD)
        print(
            "Filtering input MT by depth: DP="
            + str(min_dp)
            + ", genotype quality: GQ="
            + str(min_gq)
            + ", VAF: VAF="
            + str(min_vaf)
        )
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
    mt_with_sampleqc.write(mt_qc_outfile, overwrite=True)
    mt_with_sampleqc.cols().write(ht_qc_cols_outfile, overwrite=True)
    # stratify by pop
    pop_ht = hl.read_table(ht_qc_cols_outfile)
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

    pop_ht = pop_ht.annotate_globals(**pop_filter_ht.index_globals())
    pop_ht = pop_ht.annotate(**pop_filter_ht[pop_ht.key]).persist()
    checkpoint = pop_ht.aggregate(hl.agg.count_where(hl.len(pop_ht.qc_metrics_filters) == 0))

    print(f"=== {checkpoint} of {n_samples} exome samples found passing pop filtering")
    pop_ht.write(qc_filter_file, overwrite=True)
    pop_ht.export(sampleqc_report, delimiter="\t")

    output_globals_json = sampleqc_report + ".globals.json"
    pop_ht.globals.export(output_globals_json)


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

    # annotate mt with runid and pop
    raw_mt_file = os.path.join(mtdir, "gatk_unprocessed.mt")
    pop_ht_file = os.path.join(mtdir, "pop_assignments.ht")
    annotated_mt_file = os.path.join(mtdir, "gatk_unprocessed_with_pop.mt")
    annotate_mt("file://" + raw_mt_file, "file://" + pop_ht_file, "file://" + annotated_mt_file)

    # run sample QC and stratify by population
    mt_qc_outfile = os.path.join(mtdir, "mt_pops_sampleqc.mt")
    ht_qc_cols_outfile = os.path.join(mtdir, "mt_pops_sampleqc.ht")
    qc_filter_file = os.path.join(mtdir, "mt_pops_QC_filters.ht")
    sampleqc_report = os.path.join(annotdir, "sample_qc_by_pop.tsv")
    stratified_sample_qc(
        "file://" + annotated_mt_file,
        "file://" + mt_qc_outfile,
        "file://" + ht_qc_cols_outfile,
        "file://" + qc_filter_file,
        "file://" + sampleqc_report,
    )
    hail_utils.stop_hl(sc)


if __name__ == "__main__":
    main()
