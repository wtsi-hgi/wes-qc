"""
This is the alternative script to make population PAC
new approach that makes pruning after all filtering.

It doesn't use identification and pruning of the related samples.

HOW TO:
0. Run the script 1-hard_filter_sex_annotation
1. Run this script. You don't need to run the 2-prune_related_samples.py script.
   In will generate the
2. Run the script 3-population_pca_prediction --pca --pca-plot --assign_pops --pca-plot-assigned


Script input:
Input matrixtable (mt) is "mt_sex_annotated.mt"
Input 1000G matrix (kg_mt) is "kg_wes_regions.mt"
long_range_ld_file is "long_ld_regions.hg38.bed"
pops_file is "igsr_samples.tsv"

Script output:
"merged_with_1kg.mt" - checkpoint for the merged matrix
"merged_with_1kg_filtered.mt" - final output of the script
The matrixtable name is not 100% correct.
It was chosen to keep unchanged the process in script 3
"""

import os
from typing import Any

import hail as hl
import argparse
from gnomad.sample_qc.ancestry import assign_population_pcs, pc_project
from utils.utils import parse_config
from wes_qc import hail_utils, filtering


def merge_and_prune(
    mt: hl.MatrixTable,
    kg_mt: hl.MatrixTable,
    long_range_ld_file: str,
    pops_file: str,
    merged_mt_file: str,
    pruned_mt_file: str,
) -> None:
    # filter both matrices
    mt_filtered = filtering.filter_matrix_for_ldprune(mt, long_range_ld_file)
    kg_filtered = filtering.filter_matrix_for_ldprune(kg_mt, long_range_ld_file)
    # removing and adding needed entries to replicate filtered_mt_file structure
    mt_filtered = mt_filtered.drop(
        "AD", "DP", "GQ", "MIN_DP", "PGT", "PID", "PL", "PS", "SB", "RGQ", "callrate", "f_stat", "is_female"
    )
    kg_filtered = kg_filtered.select_entries(kg_filtered.GT)
    # merging matrices
    mt_merged = mt_filtered.union_cols(kg_filtered)
    # saving the merged matrix
    mt_merged = mt_merged.checkpoint(merged_mt_file, overwrite=True)

    # prunning of the linked Variants
    pruned_ht = hl.ld_prune(mt_merged.GT, r2=0.2)
    pruned_mt = mt_merged.filter_rows(hl.is_defined(pruned_ht[mt_merged.row_key]))
    pruned_mt = pruned_mt.select_entries(GT=hl.unphased_diploid_gt_index_call(pruned_mt.GT.n_alt_alleles()))
    # annotating known poppulations
    cohorts_pop = hl.import_table(pops_file, delimiter="\t").key_by("Sample name")
    pruned_mt = pruned_mt.annotate_cols(known_pop=cohorts_pop[pruned_mt.s]["Superpopulation code"])
    # saving matrix
    pruned_mt.write(pruned_mt_file, overwrite=True)


def main() -> None:
    # set up
    inputs = parse_config()

    data_root: str = inputs["data_root"]
    dataset_name: str = inputs["dataset_name"]
    mtdir: str = os.path.join(data_root, inputs["matrixtables_lustre_dir"])
    resourcedir = os.path.join(data_root, inputs["resource_dir"])
    annot_dir: str = os.path.join(data_root, inputs["annotation_lustre_dir"])
    plotdir = os.path.join(data_root, inputs["plots_lustre_dir"])

    # initialise hail
    sc = hail_utils.init_hl(inputs["tmp_dir"])

    mtfile = os.path.join(mtdir, "mt_sex_annotated.mt")
    mt = hl.read_matrix_table("file://" + mtfile)
    kg_mt_file = os.path.join(mtdir, "kg_wes_regions.mt")
    kg_mt = hl.read_matrix_table("file://" + kg_mt_file)

    pops_file = os.path.join(resourcedir, inputs["kg_pop"])
    long_range_ld_file = os.path.join(resourcedir, "long_ld_regions.hg38.bed")
    merged_mt_file = os.path.join(mtdir, "merged_with_1kg.mt")  # intermediate checkpoint
    # The resulting file we will use for PCA
    pruned_mt_file = os.path.join(mtdir, "merged_with_1kg_filtered.mt")

    merge_and_prune(
        mt,
        kg_mt,
        "file://" + long_range_ld_file,
        "file://" + pops_file,
        "file://" + merged_mt_file,
        "file://" + pruned_mt_file,
    )

    hail_utils.stop_hl(sc)


if __name__ == "__main__":
    main()
