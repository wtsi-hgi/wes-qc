# population prediction with PCA
import os
from typing import Any

import hail as hl
import argparse
from gnomad.sample_qc.ancestry import assign_population_pcs, pc_project
from utils.utils import parse_config
from wes_qc import hail_utils, visualize


def get_options() -> Any:
    """
    Get options from the command line
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-k", "--kg-to-mt", help="convert 1kg data to matrixtable", action="store_true")
    parser.add_argument("-m", "--merge", help="merge alspac mt with 1kg mt", action="store_true")
    parser.add_argument("-f", "--filter", help="annotate and filter merged mt", action="store_true")
    parser.add_argument("-p", "--pca", help="run pca", action="store_true")
    parser.add_argument("--pca-plot", help="Plot PCA for 1000genomes", action="store_true")
    parser.add_argument("-a", "--assign-pops", help="assign populations", action="store_true")
    parser.add_argument("--pca-plot-assigned", help="Plot PCA for assigned populations", action="store_true")
    parser.add_argument("-r", "--run", help="run all steps except kg_to_mt", action="store_true")
    args = parser.parse_args()
    return args


def merge_with_1kg(pruned_mt_file: str, kg_mt_file: str, merged_mt_file: str) -> None:
    """
    merge the birth cohort WES ld-pruned data with 1kg data of known population
    :param str pruned_mt_file: ld pruned MT file
    :param str mtdir: resources directory
    :param str merged_mt_file: merged output MT file
    """
    print("Merging with 1kg data")
    mt = hl.read_matrix_table(pruned_mt_file)
    kg_mt = hl.read_matrix_table(kg_mt_file)
    print(f"=== Variations in cohort: {mt.count_rows()}")
    print(f"=== Variations in 1kg: {kg_mt.count_rows()}")
    # in order to create a union dataset the entries and column fields must be
    # the same in each dataset. The following 2 lines take care of this.
    kg_mt = kg_mt.select_entries(kg_mt.GT)
    mt = mt.drop("callrate", "f_stat", "is_female")
    # union cols gives all samples and the rows which are found in both
    mt_joined = mt.union_cols(kg_mt)
    mt_joined.write(merged_mt_file, overwrite=True)
    print(f"=== Variations in merged cohort & 1000genomes: {mt_joined.count_rows()}")


def annotate_and_filter(merged_mt_file: str, filtered_mt_file: str, pops_file: str, long_range_ld_file: str) -> None:
    """
    Annotate with known pops for 1kg samples and filter to remove long range LD regions,
    rare variants, palidromic variants, low call rate and HWE filtering
    :param str merged_mt_file: merged birth cohort wes and 1kg MT file
    :param str filtered_mt_file: merged birth cohort wes and 1kg MT file annotated with pops and filtered
    """
    print("Adding population annotation for 1kg samples")
    mt = hl.read_matrix_table(merged_mt_file)

    # The following is 1kg superpop
    cohorts_pop = hl.import_table(pops_file, delimiter="\t").key_by("Sample name")
    mt = mt.annotate_cols(known_pop=cohorts_pop[mt.s]["Superpopulation code"])

    print("Filtering variants")
    mt_vqc = hl.variant_qc(mt, name="variant_QC_Hail")
    mt_vqc_filtered = mt_vqc.filter_rows(
        (mt_vqc.variant_QC_Hail.call_rate >= 0.99)
        & (mt_vqc.variant_QC_Hail.AF[1] >= 0.05)
        & (mt_vqc.variant_QC_Hail.p_value_hwe >= 1e-5)
    )

    long_range_ld_to_exclude = hl.import_bed(long_range_ld_file, reference_genome="GRCh38")
    mt_vqc_filtered = mt_vqc_filtered.filter_rows(
        hl.is_defined(long_range_ld_to_exclude[mt_vqc_filtered.locus]), keep=False
    )
    mt_non_pal = mt_vqc_filtered.filter_rows(
        (mt_vqc_filtered.alleles[0] == "G") & (mt_vqc_filtered.alleles[1] == "C"), keep=False
    )
    mt_non_pal = mt_non_pal.filter_rows((mt_non_pal.alleles[0] == "C") & (mt_non_pal.alleles[1] == "G"), keep=False)
    mt_non_pal = mt_non_pal.filter_rows((mt_non_pal.alleles[0] == "A") & (mt_non_pal.alleles[1] == "T"), keep=False)
    mt_non_pal = mt_non_pal.filter_rows((mt_non_pal.alleles[0] == "T") & (mt_non_pal.alleles[1] == "A"), keep=False)

    mt_non_pal.write(filtered_mt_file, overwrite=True)
    print(f"=== Variations after all filtrations : {mt_non_pal.count_rows()}")


def run_pca(
    filtered_mt_file: str,
    pca_scores_file: str,
    pca_1kg_scores_file: str,
    pca_1kg_loadings_file: str,
    pca_1kg_evals_file: str,
) -> None:
    """
    Run PCA before population prediction
    :param str filtered_mt_file: merged birth cohort wes and 1kg MT file annotated with pops and filtered
    :param str pca_sores_file: PCA scores HT file
    :param str pca_loadings_file: PCA scores HT file
    :param str pca_evals_file: PCA scores HT file
    """
    mt = hl.read_matrix_table(filtered_mt_file)
    # divide matrix to make a projection
    mt_kg = mt.filter_cols(hl.is_defined(mt.known_pop))
    mt_study = mt.filter_cols(hl.is_missing(mt.known_pop))
    # PCA for 1000 Genomes
    pca_evals, pca_scores, pca_loadings = hl.hwe_normalized_pca(mt_kg.GT, k=10, compute_loadings=True)
    pca_scores = pca_scores.annotate(known_pop=mt_kg.cols()[pca_scores.s].known_pop)
    pca_af_ht = mt_kg.annotate_rows(pca_af=hl.agg.mean(mt_kg.GT.n_alt_alleles()) / 2).rows()
    pca_loadings = pca_loadings.annotate(pca_af=pca_af_ht[pca_loadings.key].pca_af)
    # saving files
    pca_scores.write(pca_1kg_scores_file, overwrite=True)
    pca_loadings.write(pca_1kg_loadings_file, overwrite=True)
    with open(pca_1kg_evals_file, "w") as f:
        for val in pca_evals:
            f.write(str(val) + "\n")
    pca_scores = pca_scores.drop(pca_scores.known_pop)
    # projection of samples on precomputed PCs and combining of two PCA_scores tables
    print("=== Projecting PCA scores to the cohort samples ===")
    projection_PCA_scores = pc_project(mt_study, pca_loadings, loading_location="loadings", af_location="pca_af")
    union_PCA_scores = pca_scores.union(projection_PCA_scores)
    union_PCA_scores = union_PCA_scores.annotate(known_pop=mt.cols()[union_PCA_scores.s].known_pop)
    # saving union_PCA_scores
    union_PCA_scores.write(pca_scores_file, overwrite=True)


def predict_pops(pca_scores_file: str, pop_ht_file: str) -> None:
    """
    Predict populations from PCA scores
    :param str pca_sores_file: PCA scores HT file
    :param str pop_ht_file: predicted population HT file
    """
    pca_scores = hl.read_table(pca_scores_file)
    known_col = "known_pop"
    print("=== Predicting populations from PCA scores ===")
    pop_ht, pop_clf = assign_population_pcs(
        pca_scores, pca_scores.scores, known_col=known_col, n_estimators=100, prop_train=0.8, min_prob=0.5
    )
    pop_ht.write(pop_ht_file, overwrite=True)
    # This is a work around for hail <0.2.88 - convert hail table to pandas df then run assign_population_pcs
    # pop_pca_scores = pca_scores.select(known_col, pca_scores=pca_scores.scores)
    # pop_pc_pd = pop_pca_scores.to_pandas()
    # pop_pc_pd = expand_pd_array_col(pop_pc_pd, "pca_scores", 10, 'PC')

    # pop_pca_scores = pca_scores.select(known_col, pca_scores=pc_cols)
    # pop_pc_pd = pop_pca_scores.to_pandas()
    # pc_cols = [f"PC{i+1}" for i in range(10)]
    # pop_pd, pop_clf = assign_population_pcs(pop_pc_pd, pc_cols, known_col=known_col, n_estimators=100, prop_train=0.8, min_prob=0.5)


def main() -> None:
    # set up
    args = get_options()
    inputs = parse_config()

    data_root: str = inputs["data_root"]
    dataset_name: str = inputs["dataset_name"]
    mtdir: str = os.path.join(data_root, inputs["matrixtables_lustre_dir"])
    resourcedir = os.path.join(data_root, inputs["resource_dir"])
    annot_dir: str = os.path.join(data_root, inputs["annotation_lustre_dir"])
    plotdir = os.path.join(data_root, inputs["plots_lustre_dir"])

    n_pca = int(inputs["n_pca"])

    kg_mt_file = os.path.join(mtdir, "kg_wes_regions.mt")

    # initialise hail
    sc = hail_utils.init_hl(inputs["tmp_dir"])

    # combine with 1KG data
    pruned_mt_file = os.path.join(mtdir, "mt_ldpruned.mt")
    merged_mt_file = os.path.join(mtdir, "merged_with_1kg.mt")
    if args.merge or args.run:
        merge_with_1kg("file://" + pruned_mt_file, "file://" + kg_mt_file, "file://" + merged_mt_file)

    # annotate and filter
    filtered_mt_file = os.path.join(mtdir, "merged_with_1kg_filtered.mt")
    pops_file = os.path.join(resourcedir, inputs["kg_pop"])
    long_range_ld_file = os.path.join(resourcedir, "long_ld_regions.hg38.bed")
    if args.filter or args.run:
        annotate_and_filter(
            "file://" + merged_mt_file,
            "file://" + filtered_mt_file,
            "file://" + pops_file,
            "file://" + long_range_ld_file,
        )

    # run pca
    pca_scores_file = os.path.join(mtdir, "pca_scores_after_pruning.ht")
    pca_1kg_scores_file = os.path.join(mtdir, "pca_scores_1kg.ht")
    pca_1kg_loadings_file = os.path.join(mtdir, "pca_loadings_1kg.ht")
    pca_1kg_evals_file = os.path.join(annot_dir, "pca_evals_1kg.txt")  # text file may need to be without file///
    if args.pca or args.run:
        run_pca(
            "file://" + filtered_mt_file,
            "file://" + pca_scores_file,
            "file://" + pca_1kg_scores_file,
            "file://" + pca_1kg_loadings_file,
            pca_1kg_evals_file,
        )

    if args.pca_plot or args.run:
        print(f"Plotting PCA components for {pca_scores_file}")
        pca_scores = hl.read_table("file://" + pca_scores_file)
        visualize.plot_pca_bokeh(
            pca_scores, os.path.join(plotdir, f"PCA_scores_after_pruning_PCA{n_pca}.html"), n_pca, pop="known_pop"
        )

        print(f"Plotting PCA components for {pca_1kg_scores_file}")
        pca_1kg_scores = hl.read_table("file://" + pca_1kg_scores_file)
        cohorts_pop = hl.import_table("file://" + pops_file, delimiter="\t").key_by("Sample name")
        pca_1kg_scores = pca_1kg_scores.annotate(known_pop=cohorts_pop[pca_1kg_scores.s]["Superpopulation code"])
        visualize.plot_pca_bokeh(
            pca_1kg_scores, os.path.join(plotdir, f"PCA_1kg_scores_PCA{n_pca}.html"), n_pca, pop="known_pop"
        )

    # assign pops
    pop_ht_file = os.path.join(mtdir, "pop_assignments.ht")
    if args.assign_pops or args.run:
        predict_pops("file://" + pca_scores_file, "file://" + pop_ht_file)

    if args.pca_plot_assigned or args.run:
        print(f"Plotting PCA components for assigned populations: {pca_scores_file}")
        pop_ht = hl.read_table("file://" + pop_ht_file)
        pop_ht = pop_ht.transmute(scores=pop_ht.pca_scores)
        visualize.plot_pca_bokeh(
            pop_ht, os.path.join(plotdir, f"PCA_assigned_populations_{n_pca}.html"), n_pca, pop="pop"
        )

    hail_utils.stop_hl(sc)


if __name__ == "__main__":
    main()
