"""
This is the new approach that makes pruning after all filtering

HOWTO:
1. Run this script instead of script 2-prune_related_samples
2. Then, run script
3-population_pca_prediction --pca --pca-plot --assign_pops --pca-plot-assigned



In the main you need to run merge_and_prune.
merge_and_prune will use filter_matrix to filter both matrices

Input matrixtable (mt) is "mt_sex_annotated.mt"
Input 1000G matrix (kg_mt) is "kg_wes_regions.mt"
long_range_ld_file is "long_ld_regions.hg38.bed"
pops_file is "igsr_samples.tsv"
merged_mt_file name of a file for the merged matrix checkpoint
pruned_mt_file name of a file to save the pruned matrix
Rsult: pruned_mt_file is a replacement for filtered_mt_file in the PCA function in the 3rd script.


"""

import os
from typing import Any

import hail as hl
import argparse
from gnomad.sample_qc.ancestry import assign_population_pcs, pc_project
from utils.utils import parse_config
from wes_qc import hail_utils, visualize


def filter_matrix(mt: hl.MatrixTable, long_range_ld_file: str) -> hl.MatrixTable:
    # use only autosomes
    mt = mt.filter_rows(mt.locus.in_autosome())
    # split multiallelic variants and remove them
    mt = hl.split_multi_hts(mt)  # this shouldn't do anything as only biallelic sites are used
    mt = mt.filter_rows(mt.was_split == True, keep=False)
    # keep only SNPs
    mt = mt.filter_rows(hl.is_snp(mt.alleles[0], mt.alleles[1]))
    # keep good variants using hail variant_qc and thre filters
    mt_vqc = hl.variant_qc(mt, name="variant_QC_Hail")
    mt_vqc_filtered = mt_vqc.filter_rows(
        (mt_vqc.variant_QC_Hail.call_rate >= 0.99)
        & (mt_vqc.variant_QC_Hail.AF[1] >= 0.05)
        & (mt_vqc.variant_QC_Hail.p_value_hwe >= 1e-5)
    )
    # remove long ld regions
    long_range_ld_to_exclude = hl.import_bed(long_range_ld_file, reference_genome="GRCh38")
    mt_vqc_filtered = mt_vqc_filtered.filter_rows(
        hl.is_defined(long_range_ld_to_exclude[mt_vqc_filtered.locus]), keep=False
    )
    # remove palindromes
    mt_non_pal = mt_vqc_filtered.filter_rows(
        (mt_vqc_filtered.alleles[0] == "G") & (mt_vqc_filtered.alleles[1] == "C"), keep=False
    )
    mt_non_pal = mt_non_pal.filter_rows((mt_non_pal.alleles[0] == "C") & (mt_non_pal.alleles[1] == "G"), keep=False)
    mt_non_pal = mt_non_pal.filter_rows((mt_non_pal.alleles[0] == "A") & (mt_non_pal.alleles[1] == "T"), keep=False)
    mt_non_pal = mt_non_pal.filter_rows((mt_non_pal.alleles[0] == "T") & (mt_non_pal.alleles[1] == "A"), keep=False)

    return mt_non_pal


def merge_and_prune(
    mt: hl.MatrixTable,
    kg_mt: hl.MatrixTable,
    long_range_ld_file: str,
    pops_file: str,
    merged_mt_file: str,
    pruned_mt_file: str,
) -> None:
    # filter both matrices
    mt_filtered = filter_matrix(mt, long_range_ld_file)
    kg_filtered = filter_matrix(kg_mt, long_range_ld_file)
    # removing and adding needed entries to replicate filtered_mt_file structure
    mt_filtered = mt_filtered.drop(
        "AD", "DP", "GQ", "MIN_DP", "PGT", "PID", "PL", "PS", "SB", "RGQ", "callrate", "f_stat", "is_female"
    )
    kg_filtered = kg_filtered.select_entries(kg_filtered.GT)
    # merging matrices
    mt_merged = mt_filtered.union_cols(kg_filtered)
    # saving the merged matrix
    mt_merged = mt_merged.checkpoint(merged_mt_file, overwrite=True)
    # prunning
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
