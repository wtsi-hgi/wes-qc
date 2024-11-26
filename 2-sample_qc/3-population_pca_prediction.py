# population prediction with PCA
import argparse
from typing import Tuple

import hail as hl
import pandas as pd
from gnomad.sample_qc.ancestry import assign_population_pcs, pc_project

from utils.config import path_local, path_spark
from utils.utils import parse_config
from wes_qc import hail_utils, filtering, visualize


def merge_1kg_and_ldprune(
    mt: hl.MatrixTable,
    kg_mt: hl.MatrixTable,
    long_range_ld_file: str,
    merged_filtered_mt_outfile: str,
    r2_threshold: float,
    call_rate_threshold: float,
    af_threshold: float,
    hwe_threshold: float,
    **kwargs,
) -> hl.MatrixTable:
    """
    Merge input and 1kg matrix and run LD pruning
    :param mt: input matrix table
    :param kg_mt: 1kg matrix table
    :param long_range_ld_file: Long range LD file
    :param merged_filtered_mt_outfile: Merged and filtered matrix output file
    :param r2_threshold: Correlation threshold for LD pruning
    :param call_rate_threshold: Call rate threshold for filtering
    :param af_threshold: Allele frequency threshold for filtering
    :param hwe_threshold: Hardy-Weinberg equilibrium threshold for filtering
    """

    # filter sample matrix
    mt_filtered = filtering.filter_matrix_for_ldprune(
        mt, path_spark(long_range_ld_file), call_rate_threshold, af_threshold, hwe_threshold
    )

    # removing and adding needed entries to replicate filtered_mt_file structure
    mt_filtered = mt_filtered.drop(
        "AD", "DP", "GQ", "MIN_DP", "PGT", "PID", "PL", "PS", "SB", "RGQ", "callrate", "f_stat", "is_female"
    )
    mt_filtered = mt_filtered.annotate_cols(known_pop=hl.null(hl.tstr))
    # merging matrices
    mt_merged = mt_filtered.union_cols(kg_mt)
    # saving the merged matrix
    mt_merged = mt_merged.checkpoint(merged_filtered_mt_outfile, overwrite=True)

    # prunning of the linked Variants
    pruned_ht = hl.ld_prune(mt_merged.GT, r2=r2_threshold)
    pruned_mt = mt_merged.filter_rows(hl.is_defined(pruned_ht[mt_merged.row_key]))
    pruned_mt = pruned_mt.select_entries(GT=hl.unphased_diploid_gt_index_call(pruned_mt.GT.n_alt_alleles()))

    return pruned_mt


def pop_pca(mt: hl.MatrixTable, pca_components: int, pca_1kg_evals_file: str, **kwargs) -> Tuple[hl.Table, hl.Table]:
    """
    Run PCA before population prediction
    :param mt: merged birth cohort wes and 1kg MT file annotated with pops and filtered
    :param str pca_evals_outfile: PCA scores HT file
    :param int pca_components: the number of principal components
    """

    # divide matrix to make a projection
    # TODO: NEED to fix joining step to make his work
    mt_kg = mt.filter_cols(hl.is_defined(mt.known_pop))  # The golden source 1000G with known population
    mt_study = mt.filter_cols(hl.is_missing(mt.known_pop))
    # PCA for golden source 1000 Genomes
    pca_1kg_evals, pca_1kg_scores, pca_1kg_loadings = hl.hwe_normalized_pca(
        mt_kg.GT, k=pca_components, compute_loadings=True
    )
    pca_1kg_scores = pca_1kg_scores.annotate(known_pop=mt_kg.cols()[pca_1kg_scores.s].known_pop)
    pca_af_ht = mt_kg.annotate_rows(pca_af=hl.agg.mean(mt_kg.GT.n_alt_alleles()) / 2).rows()
    pca_1kg_loadings = pca_1kg_loadings.annotate(pca_af=pca_af_ht[pca_1kg_loadings.key].pca_af)

    with open(pca_1kg_evals_file, "w") as f:
        for val in pca_1kg_evals:
            f.write(str(val) + "\n")

    pca_1kg_scores_nopop = pca_1kg_scores.drop(pca_1kg_scores.known_pop)

    # projection of samples on precomputed PCs and combining of two PCA_scores tables
    print("=== Projecting PCA scores to the cohort samples ===")
    projection_PCA_scores = pc_project(mt_study, pca_1kg_loadings, loading_location="loadings", af_location="pca_af")
    union_PCA_scores = pca_1kg_scores_nopop.union(projection_PCA_scores)
    union_PCA_scores = union_PCA_scores.annotate(known_pop=mt.cols()[union_PCA_scores.s].known_pop)

    return pca_1kg_scores, pca_1kg_loadings, union_PCA_scores


def append_row(df, row):
    return pd.concat([df, pd.DataFrame([row], columns=row.index)]).reset_index(drop=True)


def predict_pops(pop_pca_scores, gnomad_pc_n_estimators, gnomad_prop_train, gnomad_min_prob, pop_ht_outtsv, **kwargs):
    """
    Predict populations from PCA scores
    :param str pca_sores_file: PCA scores HT file
    :param str pop_ht_file: predicted population HT file
    :param str pop_ht_tsv: population tsv file
    """
    known_col = "known_pop"
    pop_ht, pop_clf = assign_population_pcs(
        pop_pca_scores,
        pop_pca_scores.scores,
        known_col=known_col,
        n_estimators=gnomad_pc_n_estimators,
        prop_train=gnomad_prop_train,
        min_prob=gnomad_min_prob,
    )

    # convert to pandas and put in only pops files, add excluded sample back
    pop_ht_df = pop_ht.to_pandas()
    pop_ht_df2 = pop_ht_df[["s", "pop"]]
    # new_row = pd.Series({"s": "EGAN00004311029", "pop": "oth"}) # TODO: Why we're including hardcoded string. Commented for now
    # pop_ht_df2 = append_row(pop_ht_df2, new_row)

    pop_ht_tsv = path_local(pop_ht_outtsv)
    pop_ht_df2.to_csv(pop_ht_tsv, sep="\t")
    return pop_ht


def get_options() -> argparse.Namespace:
    """
    Get options from the command line
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--merge-and-ldprune", help="Merge with 1000G data and run LD-pruning", action="store_true")
    parser.add_argument("--pca", help="Run pca", action="store_true")
    parser.add_argument("--pca-plot", help="Plot PCA for 1000genomes", action="store_true")
    parser.add_argument("--assign_pops", help="Assign populations", action="store_true")
    parser.add_argument("--pca-plot-assigned", help="Plot PCA for assigned populations", action="store_true")
    parser.add_argument("--all", help="Run all steps", action="store_true")
    args = parser.parse_args()
    return args


def main():
    # = STEP SETUP = #
    config = parse_config()
    args = get_options()
    if args.all:
        args.merge_and_ldprune = True
        args.pca = True
        args.pca_plot = True
        args.assign_pops = True
        args.plot_pca_assigned = True

    tmp_dir = config["general"]["tmp_dir"]

    # = STEP PARAMETERS = #
    ## No parameters for this step

    # = STEP DEPENDENCIES = #
    mtfile = path_spark(config["step2"]["impute_sex"]["sex_mt_outfile"])
    kg_mt_file = path_spark(config["step1"]["create_1kg_mt"]["kg_out_mt"])
    kg_pop_file = path_spark(config["step1"]["create_1kg_mt"]["kg_pop_file"])

    # = STEP OUTPUTS = #
    pruned_mt_file = path_spark(config["step2"]["merge_1kg_and_ldprune"]["filtered_and_pruned_mt_outfile"])
    pca_1kg_scores_file = path_spark(config["step2"]["pop_pca"]["pca_1kg_scores_file"])
    pca_1kg_loadings_file = path_spark(config["step2"]["pop_pca"]["pca_1kg_loadings_file"])
    pca_union_scores_file = path_spark(config["step2"]["pop_pca"]["pca_union_scores_file"])

    # = STEP LOGIC = #
    _ = hail_utils.init_hl(tmp_dir)

    if args.merge_and_ldprune:
        mt = hl.read_matrix_table(mtfile)
        kg_mt = hl.read_matrix_table(kg_mt_file)
        pruned_mt = merge_1kg_and_ldprune(mt, kg_mt, **config["step2"]["merge_1kg_and_ldprune"])
        pruned_mt.write(pruned_mt_file, overwrite=True)

    # run pca
    if args.pca:
        filtered_mt = hl.read_matrix_table(pruned_mt_file)
        pca_1kg_scores, pca_1kg_loadings, union_PCA_scores = pop_pca(filtered_mt, **config["step2"]["pop_pca"])
        pca_1kg_scores.write(pca_1kg_scores_file, overwrite=True)
        pca_1kg_loadings.write(pca_1kg_loadings_file, overwrite=True)
        union_PCA_scores.write(pca_union_scores_file, overwrite=True)

    if args.pca_plot:
        print(f"Plotting PCA components for {pca_union_scores_file}")
        pca_union_scores = hl.read_table(pca_union_scores_file)
        print("Union components:", pca_union_scores.count())
        pop_pca_union_graph = config["step2"]["plot_pop_pca"]["pop_pca_union_graph"]
        n_pca = config["step2"]["plot_pop_pca"]["pca_components"]
        visualize.plot_pop_pca(pca_union_scores, pop_pca_union_graph, n_pca, pop="known_pop")

        print(f"Plotting PCA components for {pca_1kg_scores_file}")
        pop_pca_1kg_graph = config["step2"]["plot_pop_pca"]["pop_pca_1kg_graph"]
        pca_1kg_scores = hl.read_table(pca_1kg_scores_file)
        print("1kg components:", pca_1kg_scores.count())
        cohorts_pop = hl.import_table(kg_pop_file, delimiter="\t").key_by("Sample name")
        pca_1kg_scores = pca_1kg_scores.annotate(known_pop=cohorts_pop[pca_1kg_scores.s]["Superpopulation code"])
        visualize.plot_pop_pca(pca_1kg_scores, pop_pca_1kg_graph, n_pca, pop="known_pop")

    # assign pops
    pop_ht_file = path_spark(config["step2"]["predict_pops"]["pop_ht_outfile"])
    if args.assign_pops:
        union_PCA_scores = hl.read_table(pca_union_scores_file)
        pop_ht = predict_pops(union_PCA_scores, **config["step2"]["predict_pops"])
        pop_ht.write(pop_ht_file, overwrite=True)

    if args.pca_plot_assigned:
        n_pca = config["step2"]["plot_pop_pca_assigned"]["pca_components"]

        print(f"Plotting PCA components for assigned populations: {pop_ht_file}")
        pop_ht = hl.read_table(pop_ht_file)
        pop_ht = pop_ht.transmute(scores=pop_ht.pca_scores)
        print(f"Total samples: {pop_ht.count()}")
        pop_assigned_pca_union_graph = config["step2"]["plot_pop_pca_assigned"]["pop_assigned_pca_union_graph"]

        visualize.plot_pop_pca(pop_ht, pop_assigned_pca_union_graph, n_pca, pop="pop")

        print("Plotting PCA components for the dataset:")
        pop_ht_datasetonly = pop_ht.filter(~hl.is_defined(pop_ht.known_pop))
        print(f"Total samples: {pop_ht_datasetonly.count()}")
        pop_pca_assigned_1kg_graph = config["step2"]["plot_pop_pca_assigned"]["pop_assigned_pca_dataset_graph"]
        visualize.plot_pop_pca(
            pop_ht_datasetonly,
            pop_pca_assigned_1kg_graph,
            n_pca,
            pop="pop",
        )
