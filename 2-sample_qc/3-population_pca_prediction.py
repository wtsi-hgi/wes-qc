# population prediction with PCA
import hail as hl
import argparse
import pandas as pd
from gnomad.sample_qc.ancestry import assign_population_pcs
from utils.utils import parse_config
from utils.config import path_local, path_spark
from wes_qc import hail_utils, filtering


def merge_1kg_and_ldprune(mt: hl.MatrixTable, kg_mt: hl.MatrixTable, conf: dict) -> hl.MatrixTable:
    # filter sample matrix
    long_range_ld_file = conf["long_range_ld_file"]
    merged_1kg_filterd_mt_file = conf["merged_filtered_mt_outfile"]
    mt_filtered = filtering.filter_matrix_for_ldprune(mt, long_range_ld_file)
    # TODO: probably next step is obsolete because we filtered 1kg matrix in the step 1.4
    kg_filtered = filtering.filter_matrix_for_ldprune(kg_mt, long_range_ld_file)

    # removing and adding needed entries to replicate filtered_mt_file structure
    mt_filtered = mt_filtered.drop(
        "AD", "DP", "GQ", "MIN_DP", "PGT", "PID", "PL", "PS", "SB", "RGQ", "callrate", "f_stat", "is_female"
    )
    kg_filtered = kg_filtered.select_entries(kg_filtered.GT)
    # merging matrices
    mt_merged = mt_filtered.union_cols(kg_filtered)
    # saving the merged matrix
    mt_merged = mt_merged.checkpoint(merged_1kg_filterd_mt_file, overwrite=True)

    # prunning of the linked Variants
    pruned_ht = hl.ld_prune(mt_merged.GT, r2=0.2)
    pruned_mt = mt_merged.filter_rows(hl.is_defined(pruned_ht[mt_merged.row_key]))
    pruned_mt = pruned_mt.select_entries(GT=hl.unphased_diploid_gt_index_call(pruned_mt.GT.n_alt_alleles()))

    # saving matrix
    return pruned_mt


def run_pca(filtered_mt_file: hl.MatrixTable, config: dict):
    """
    Run PCA before population prediction
    :param str filtered_mt_file: merged birth cohort wes and 1kg MT file annotated with pops and filtered
    :param str pca_sores_file: PCA scores HT file
    :param str pca_loadings_file: PCA scores HT file
    :param str pca_evals_file: PCA scores HT file
    """
    conf = config["step2"]["run_pca"]
    # mt = hl.read_matrix_table(filtered_mt_file)
    mt = filtered_mt_file  # FIXME: testing accepting mt
    # exclude EGAN00004311029 this is a huge outlier on PCA and skews all the PCs
    to_exclude = ["EGAN00004311029"]  # TODO: add as a list parameter to the config
    mt = mt.filter_cols(~hl.set(to_exclude).contains(mt.s))

    # TODO: move k to config
    pca_evals, pca_scores, pca_loadings = hl.hwe_normalized_pca(mt.GT, k=4, compute_loadings=True)
    pca_scores = pca_scores.annotate(known_pop=mt.cols()[pca_scores.s].known_pop)
    pca_scores_file = path_spark(conf["pca_scores_outfile"])
    pca_scores.write(pca_scores_file, overwrite=True)
    pca_loadings_file = path_spark(conf["pca_loadings_outfile"])
    pca_loadings.write(pca_loadings_file, overwrite=True)
    pca_evals_file = path_local(conf["pca_evals_outfile"])

    with open(pca_evals_file, "w") as f:
        for val in pca_evals:
            f.write(str(val) + "\n")

    # don't return anything as there are multiple outputs


def append_row(df, row):
    return pd.concat([df, pd.DataFrame([row], columns=row.index)]).reset_index(drop=True)


def predict_pops(config: dict):
    """
    Predict populations from PCA scores
    :param str pca_sores_file: PCA scores HT file
    :param str pop_ht_file: predicted population HT file
    :param str pop_ht_tsv: population tsv file
    """
    conf = config["step2"]["predict_pops"]
    pca_scores_file = path_spark(conf["pca_scores_file"])
    pca_scores = hl.read_table(pca_scores_file)
    known_col = "known_pop"
    pop_ht, pop_clf = assign_population_pcs(
        pca_scores,
        pca_scores.scores,
        known_col=known_col,
        n_estimators=conf["gnomad_pc_n_estimators"],
        prop_train=conf["gnomad_prop_train"],
        min_prob=conf["gnomad_min_prob"],
    )
    pop_ht_file = path_spark(conf["pop_ht_outfile"])
    pop_ht.write(pop_ht_file, overwrite=True)
    # convert to pandas and put in only pops files, add excluded sample back
    pop_ht_df = pop_ht.to_pandas()
    pop_ht_df2 = pop_ht_df[["s", "pop"]]
    new_row = pd.Series({"s": "EGAN00004311029", "pop": "oth"})
    pop_ht_df2 = append_row(pop_ht_df2, new_row)
    print(pop_ht_df2)

    pop_ht_tsv = path_local(conf["pop_ht_outtsv"])
    print(pop_ht_tsv)

    pop_ht_df2.to_csv(pop_ht_tsv, sep="\t")

    # This is a work around for hail <0.2.88 - convert hail table to pandas df then run assign_population_pcs
    # pop_pca_scores = pca_scores.select(known_col, pca_scores=pca_scores.scores)
    # pop_pc_pd = pop_pca_scores.to_pandas()
    # pop_pc_pd = expand_pd_array_col(pop_pc_pd, "pca_scores", 10, 'PC')

    # pop_pca_scores = pca_scores.select(known_col, pca_scores=pc_cols)
    # pop_pc_pd = pop_pca_scores.to_pandas()
    # pc_cols = [f"PC{i+1}" for i in range(10)]
    # pop_pd, pop_clf = assign_population_pcs(pop_pc_pd, pc_cols, known_col=known_col, n_estimators=100, prop_train=0.8, min_prob=0.5)


def get_options() -> argparse.Namespace:
    """
    Get options from the command line
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--merge-and-ldprune", help="Merge with 1000G data and run LD-pruning", action="store_true")
    parser.add_argument("--pca", help="Run pca", action="store_true")
    parser.add_argument("--assign_pops", help="Assign populations", action="store_true")
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
        args.assign_pops = True

    tmp_dir = config["general"]["tmp_dir"]

    # = STEP PARAMETERS = #
    conf = config["step2"][""]

    # = STEP DEPENDENCIES = #

    # = STEP OUTPUTS = #

    # = STEP LOGIC = #
    _ = hail_utils.init_hl(tmp_dir)

    if args.merge_and_ldprune:
        # = STEP DEPENDENCIES = #
        mtfile = config["step2"]["impute_sex"]["sex_mt_outfile"]
        kg_mt_file = config["step1"]["create_1kg_mt"]["kg_out_mt"]
        pruned_mt_file = conf["filtered_and_pruned_mt_outfile"]

        # = STEP LOGIC = #
        _ = hail_utils.init_hl(tmp_dir)
        mt = hl.read_matrix_table("file://" + mtfile)
        kg_mt = hl.read_matrix_table("file://" + kg_mt_file)
        pruned_mt = merge_1kg_and_ldprune(mt, kg_mt, conf)
        pruned_mt.write(pruned_mt_file, overwrite=True)

    # run pca
    if args.pca or args.run:
        filtered_mt_file = path_spark(config["step2"]["merge_1kg_and_ldprune"]["filtered_and_pruned_mt_outfile"])
        filtered_mt = hl.read_matrix_table(filtered_mt_file)
        run_pca(filtered_mt, config)

    # assign pops
    if args.assign_pops or args.run:
        predict_pops(config)


if __name__ == "__main__":
    main()
