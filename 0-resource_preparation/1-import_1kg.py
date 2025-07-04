"""
This script prepares the 1000 Genome matrixtable,
starting from the published 1000G VCF files

You need to run this script only once to prepare the 1000 Genome matrix
Then you can reuse it with other datasets

To prepare data for this script
(you can use the script in scripts/1kg_download):
* Download 1000G VCFs
* Remove SV and leave only SNP and indels
"""

import argparse
from typing import Any

import hail as hl  # type: ignore

from wes_qc.hail_utils import path_spark
from wes_qc.config import parse_config
from wes_qc import filtering, hail_utils


def create_1kg_mt(vcf_indir: str, kg_pop_file: str, **kwargs: dict) -> hl.MatrixTable:
    """
    Create matrixtable of 1kg data
    :param str vcf_indir: the directory with 1KG VCF files
    :param str kg_pop_file: Assigns superpopulations
    """
    print(f"Loading VCFs from {vcf_indir}")
    objects = hl.utils.hadoop_ls(vcf_indir)
    vcfs = [vcf["path"] for vcf in objects if (vcf["path"].startswith("file") and vcf["path"].endswith("vcf.gz"))]
    # create MT
    kg_unprocessed_mt = hl.import_vcf(vcfs, array_elements_required=False, force_bgz=True)
    # Annotating known populations
    kg_pop_file = path_spark(kg_pop_file)
    cohorts_pop = hl.import_table(kg_pop_file, delimiter="\t").key_by("Sample name")
    kg_unprocessed_mt = kg_unprocessed_mt.annotate_cols(
        known_pop=cohorts_pop[kg_unprocessed_mt.s]["Superpopulation code"]
    )
    # Renaming samples to avoid name clashes with cohorts samples
    kg_unprocessed_mt = kg_unprocessed_mt.key_cols_by()
    kg_unprocessed_mt = kg_unprocessed_mt.transmute_cols(s=hl.str("1kg-for-pop-pca_") + kg_unprocessed_mt.s)
    kg_unprocessed_mt = kg_unprocessed_mt.key_cols_by(kg_unprocessed_mt.s)

    return kg_unprocessed_mt


def kg_filter_and_ldprune(
    kg_unprocessed_mt: hl.MatrixTable,
    long_range_ld_file: str,
    call_rate_threshold: float,
    af_threshold: float,
    hwe_threshold: float,
    r2_threshold: float,
    **kwargs,
) -> hl.MatrixTable:
    """
    Filter and prune the 1kg data
    :param kg_unprocessed_mt: The KG MT to filter and prune
    :param long_range_ld_file: The long range LD file
    :param call_rate_threshold: The call rate threshold
    :param af_threshold: The allele frequency threshold
    :param hwe_threshold: The HWE threshold
    :param r2_threshold: Squared correlation threshold: https://hail.is/docs/0.2/methods/genetics.html#hail.methods.ld_prune
    :return: The filtered and pruned 1kg MT
    """
    long_range_ld_file = path_spark(long_range_ld_file)
    # Filtering for good variations to make LD prune
    kg_mt_filtered = filtering.filter_matrix_for_ldprune(
        kg_unprocessed_mt,
        long_range_ld_file,
        call_rate_threshold,
        af_threshold,
        hwe_threshold,
    )
    # LD pruning - removing variation regions that are related to each other
    pruned_kg_ht = hl.ld_prune(kg_mt_filtered.GT, r2=r2_threshold)
    pruned_kg_mt = kg_mt_filtered.filter_rows(hl.is_defined(pruned_kg_ht[kg_mt_filtered.row_key]))
    pruned_kg_mt = pruned_kg_mt.select_entries(GT=hl.unphased_diploid_gt_index_call(pruned_kg_mt.GT.n_alt_alleles()))
    return pruned_kg_mt


def run_pc_relate(
    pruned_mt: hl.MatrixTable,
    relatedness_ht_file,
    scores_file: str,
    pca_components: int,
    kin_threshold: float,
    hl_pc_related_kwargs=dict(),
    **kwargs,
) -> hl.Table:
    """
    Runs PC relate on pruned MT
    :param str pruned_mt: matrixtable to prune
    :param str relatedness_ht_file: relatedness ht file
    :param str scores_file: file to wtire scores ht
    :param int pca_components:  the number of principal components
    :param dict hl_pc_related_kwargs: kwargs to pass to HL PC relate
    """
    if hl_pc_related_kwargs is None:
        hl_pc_related_kwargs = {}
    relatedness_ht_file = path_spark(relatedness_ht_file)
    scores_file = path_spark(scores_file)

    print("=== Running PC relate")
    eig, scores, _ = hl.hwe_normalized_pca(pruned_mt.GT, k=pca_components, compute_loadings=False)
    scores.write(scores_file, overwrite=True)

    print("=== Calculating relatedness (this step usually takes a while)")
    relatedness_ht = hl.pc_relate(
        pruned_mt.GT,
        scores_expr=scores[pruned_mt.col_key].scores,
        **hl_pc_related_kwargs,
    )
    relatedness_ht.write(relatedness_ht_file, overwrite=True)
    # prune individuals to be left with unrelated - creates a table containing one column - samples to remove
    pairs = relatedness_ht.filter(relatedness_ht["kin"] > kin_threshold)
    related_samples_to_remove = hl.maximal_independent_set(pairs.i, pairs.j, keep=False)
    return related_samples_to_remove


def kg_remove_related_samples(kg_mt: hl.MatrixTable, related_samples_to_remove: hl.Table) -> hl.MatrixTable:
    variants, samples = kg_mt.count()
    print(f"=== Loaded form initial table: {samples} samples, {variants} variants.")
    print("=== Removing related samples")
    kg_mt_remove_related = kg_mt.filter_cols(hl.is_defined(related_samples_to_remove[kg_mt.col_key]), keep=False)
    variants, samples = kg_mt_remove_related.count()
    print(f"=== Remains after removing related samples: {samples} samples, {variants} variants.")
    return kg_mt_remove_related


def get_options() -> Any:
    """
    Get options from the command line
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--kg-to-mt", help="Convert 1kg data to matrixtable", action="store_true")
    parser.add_argument(
        "--kg-filter-and-prune",
        help="Prune related variants from 1kg matrix",
        action="store_true",
    )
    parser.add_argument("--kg-pc-relate", help="Run PC relate for 1KG ", action="store_true")
    parser.add_argument("--kg-remove-related-samples", help="Run PC relate for 1KG", action="store_true")
    parser.add_argument("--all", help="Run All steps", action="store_true")
    args = parser.parse_args()
    return args


def main() -> None:
    # = STEP SETUP = #
    config = parse_config()
    args = get_options()
    if args.all:
        args.kg_to_mt = True
        args.kg_filter_and_prune = True
        args.kg_pc_relate = True
        args.kg_remove_related_samples = True

    tmp_dir = config["general"]["tmp_dir"]

    # = STEP PARAMETERS = #
    conf = config["step0"]["create_1kg_mt"]

    # = STEP DEPENDENCIES = #
    vcf_indir = path_spark(conf["indir"])

    # = STEP OUTPUTS = #
    kg_unprocessed_mt_file = path_spark(conf["kg_unprocessed"])
    pruned_kg_file = path_spark(conf["pruned_kg_file"])
    samples_to_remove_file = path_spark(conf["samples_to_remove_file"])
    kg_mt_file = path_spark(conf["kg_out_mt"])

    # = STEP LOGIC = #
    _ = hail_utils.init_hl(tmp_dir)

    if args.kg_to_mt:
        kg_mt = create_1kg_mt(vcf_indir, **conf)
        print(f"Saving as hail mt to {kg_unprocessed_mt_file}")
        kg_mt.write(kg_unprocessed_mt_file, overwrite=True)

    if args.kg_filter_and_prune:
        kg_unprocessed_mt = hl.read_matrix_table(kg_unprocessed_mt_file)
        pruned_kg_mt = kg_filter_and_ldprune(kg_unprocessed_mt, **conf)
        pruned_kg_mt.write(pruned_kg_file, overwrite=True)

    if args.kg_pc_relate:
        pruned_kg_mt = hl.read_matrix_table(pruned_kg_file)
        related_samples_to_remove = run_pc_relate(pruned_kg_mt, **conf)
        related_samples_to_remove.write(samples_to_remove_file, overwrite=True)

    if args.kg_remove_related_samples:
        kg_mt = hl.read_matrix_table(kg_unprocessed_mt_file)
        related_samples_to_remove = hl.read_table(samples_to_remove_file)
        kg_mt_remove_related = kg_remove_related_samples(kg_mt, related_samples_to_remove)
        kg_mt_remove_related.write(kg_mt_file, overwrite=True)


if __name__ == "__main__":
    main()
