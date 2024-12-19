# identify and prune related samples prior to PCA
# use mt with hard filters and sex annotation from 2-sample_qc/1-hard_filters_sex_annotation.py
import hail as hl
import os
from utils.utils import parse_config, path_local, path_spark
import bokeh.plotting as bkplt
import bokeh.layouts as bklayouts
from wes_qc import hail_utils


def prune_mt(mt: hl.MatrixTable, ld_prune_args, **kwargs) -> hl.MatrixTable:
    """
    Splits multiallelic sites and runs ld pruning
    Filter to autosomes before LD pruning to decrease sample size - autosomes only wanted for later steps
    :param MatrixTable mt: input MT containing variants to be pruned
    :param dict config:
    :return: Pruned MatrixTable
    :rtype: hl.MatrixTable

    `hl.ld_prune` returns a maximal subset of variants that are nearly uncorrelated within each window.
    Requires the dataset to contain only diploid genotype calls and no multiallelic variants.
    See also: https://hail.is/docs/0.2/methods/genetics.html#hail.methods.ld_prune
    """

    print("=== Filtering to autosomes")
    mt = mt.filter_rows(mt.locus.in_autosome())
    print("=== Splitting multiallelic sites")
    mt = hl.split_multi_hts(mt)  # this shouldn't do anything as only biallelic sites are used
    print("=== Performing LD pruning")
    pruned_ht = hl.ld_prune(mt.GT, **ld_prune_args)
    pruned_mt = mt.filter_rows(hl.is_defined(pruned_ht[mt.row_key]))
    pruned_mt = pruned_mt.select_entries(GT=hl.unphased_diploid_gt_index_call(pruned_mt.GT.n_alt_alleles()))
    return pruned_mt


def prune_pc_relate(
    pruned_mt: hl.MatrixTable, pca_components, pc_relate_args, relatedness_column, relatedness_threshold, **kwargs
) -> (hl.Table, hl.Table, hl.Table):
    """
    Runs PC relate on pruned MT
    :param hl.MatrixTable pruned_mt: ld pruned MT file
    :param dict config:
    :return: mt with a column of too related samples to remove
    :rtype: hl.MatrixTable

    See https://hail.is/docs/0.2/methods/relatedness.html#hail.methods.pc_relate

    ### Config fields
    step2.pc_relate.relatedness_ht_file : output path : relatedness ht file
    step2.pc_relate.samples_to_remove_file : output path : ht file with a list of samples to remove
    step2.pc_relate.scores_file : output path: : ht file with hwe_normalized_pca scores
    step2.pc_relate.pca_components : int : number of components for the hwe_normalized_pca
    step2.pc_relate.relatedness_column : str : pc_relate output column to use when applying the relatedness_threshold
    step2.pc_relate.relatedness_threshold : float : samples with relatedness_column greater than relatedness_threshold are to be removed

    ### Optional config fields
    step2.pc_relate.pc_relate_args.min_individual_maf : float
    step2.pc_relate.pc_relate_args.block_size : float
    step2.pc_relate.pc_relate_args.min_kinship : float
    step2.pc_relate.pc_relate_args.statistics : str
    step2.pc_relate.pc_relate_args.k : int
    step2.pc_relate.pc_relate_args.include_self_kinship : bool
    """

    print("=== Running PC relate")
    eig, scores, _ = hl.hwe_normalized_pca(pruned_mt.GT, k=pca_components, compute_loadings=False)
    print("=== Calculating relatedness")
    relatedness_ht = hl.pc_relate(pruned_mt.GT, scores_expr=scores[pruned_mt.col_key].scores, **pc_relate_args)
    # prune individuals to be left with unrelated - creates a table containing one column - samples to remove
    pairs = relatedness_ht.filter(relatedness_ht[relatedness_column] > relatedness_threshold)
    related_samples_to_remove = hl.maximal_independent_set(pairs.i, pairs.j, keep=False)
    return related_samples_to_remove, scores, relatedness_ht


def plot_relatedness(relatedness_ht: hl.Table, relatedness_plotfile: str, **kwargs) -> None:
    """Plot relatedness scores."""
    p1 = hl.plot.histogram(relatedness_ht.kin, title="Kinship distribution")
    p2 = hl.plot.scatter(relatedness_ht.kin, relatedness_ht.ibd2, xlabel="Kinship", ylabel="IBD2")
    layout = bklayouts.gridplot([p1, p2], ncols=2)
    bkplt.output_file(relatedness_plotfile)
    bkplt.save(layout)


# TODO: How is this step different from 2-sample_qc/3-population_pca_prediction.py/run_pca()? Only PCA plotting?
def run_population_pca(
    pruned_mt: hl.MatrixTable,
    samples_to_remove: hl.Table,
    plink_outfile,
    pca_components,
    plot_outfile,
    **kwargs,
) -> (hl.MatrixTable, hl.Table, hl.Table):
    """
    Runs PCA and creates a matrix table of non-related individuals with PCA scores
    Remove related samples from PC relate from pruned MT and run PCA
    """

    print("=== Running population PCA")
    print("=== Removing related samples")
    pca_mt = pruned_mt.filter_cols(hl.is_defined(samples_to_remove[pruned_mt.col_key]), keep=False)
    variants, samples = pca_mt.count()
    print(f"=== Survived {samples} samples after relatedness step.")

    plink_mt = pca_mt.annotate_cols(uid=pca_mt.s).key_cols_by("uid")
    hl.export_plink(dataset=plink_mt, output=path_spark(plink_outfile), fam_id=plink_mt.uid, ind_id=plink_mt.uid)

    print("=== Running PCA")
    pca_evals, pca_scores, pca_loadings = hl.hwe_normalized_pca(pca_mt.GT, k=pca_components, compute_loadings=True)
    pca_af_ht = pca_mt.annotate_rows(pca_af=hl.agg.mean(pca_mt.GT.n_alt_alleles()) / 2).rows()
    pca_loadings = pca_loadings.annotate(pca_af=pca_af_ht[pca_loadings.key].pca_af)
    pca_mt = pca_mt.annotate_cols(scores=pca_scores[pca_mt.col_key].scores)

    print("=== Plotting PC1 vs PC2")
    os.makedirs(os.path.dirname(plot_outfile), exist_ok=True)
    p = hl.plot.scatter(pca_mt.scores[0], pca_mt.scores[1], title="PCA", xlabel="PC1", ylabel="PC2")
    print(f"=== Saving Relatedness PCA plot to {plot_outfile}")
    bkplt.output_file(plot_outfile)
    bkplt.save(p)

    return pca_mt, pca_scores, pca_loadings


def main():
    # = STEP SETUP = #
    config = parse_config()
    tmp_dir = config["general"]["tmp_dir"]

    # = STEP PARAMETERS = #
    ## No parameters for this step

    # = STEP DEPENDENCIES = #
    mt_infile = config["step2"]["impute_sex"]["sex_mt_outfile"]

    # = STEP OUTPUTS = #
    mtoutfile = config["step2"]["prune"]["pruned_mt_file"]

    # = STEP LOGIC = #
    _ = hail_utils.init_hl(tmp_dir)

    # ensure plotting directory exists
    pltdir = path_local(config["general"]["plots_dir"])
    if not os.path.exists(pltdir):
        os.makedirs(pltdir)

    # load input mt
    mt = hl.read_matrix_table(path_spark(mt_infile))

    # ld prune to get a table of variants which are not correlated
    pruned_mt = prune_mt(mt, **config["step2"]["prune"])
    pruned_mt.write(path_spark(mtoutfile), overwrite=True)

    # run pcrelate
    related_samples_to_remove_ht, scores, relatedness_ht = prune_pc_relate(
        pruned_mt, **config["step2"]["prune_pc_relate"]
    )
    scores_file = config["step2"]["prune_pc_relate"]["scores_file"]
    scores.write(path_spark(scores_file), overwrite=True)  # output
    related_samples_to_remove_ht.write(
        path_spark(config["step2"]["prune_pc_relate"]["samples_to_remove_file"]), overwrite=True
    )
    related_samples_to_remove_ht.export(path_spark(config["step2"]["prune_pc_relate"]["samples_to_remove_tsv"]))
    relatedness_ht.export(path_spark(config["step2"]["prune_pc_relate"]["relatedness_outfile"]))

    # plot relatedness
    plot_relatedness(relatedness_ht, **config["step2"]["prune_pc_relate"])

    # run PCA
    pca_mt, pca_scores, pca_loadings = run_population_pca(
        pruned_mt, related_samples_to_remove_ht, **config["step2"]["prune_plot_pca"]
    )
    pca_mt.write(path_spark(config["step2"]["prune_plot_pca"]["pca_mt_file"]), overwrite=True)
    pca_scores.write(path_spark(config["step2"]["prune_plot_pca"]["pca_scores_file"]), overwrite=True)
    pca_loadings.write(path_spark(config["step2"]["prune_plot_pca"]["pca_loadings_file"]), overwrite=True)  # output


if __name__ == "__main__":
    main()
