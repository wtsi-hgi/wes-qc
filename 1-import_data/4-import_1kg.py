"""
This script prepares the 1000 Genome matrixtable,
starting from the published 1000G VCF files

You need to run this script only once to prepare the 1000 Genome matrix
Then ou can reuse it with other datasets

To prepare data for this script
(you can use the script in scripts/1kg_download):
* Download 1000G VCFs
* Remove SV and leave only SNP and indels
"""

import argparse
from typing import Any

import hail as hl  # type: ignore

from utils.utils import parse_config
from utils.config import path_spark
from wes_qc import hail_utils, filtering


def create_1kg_mt(vcf_indir: str, kg_unprocessed_mt: str) -> None:
    """
    Create matrixtable of 1kg data
    :param str resourcedir: resources directory
    :param str mtdir: matrixtable directory
    """
    print(f"Loading VCFs from {vcf_indir}")
    objects = hl.utils.hadoop_ls(vcf_indir)
    vcfs = [vcf["path"] for vcf in objects if (vcf["path"].startswith("file") and vcf["path"].endswith("vcf.gz"))]

    # create and save MT
    mt = hl.import_vcf(vcfs, array_elements_required=False, force_bgz=True)
    print(f"Saving as hail mt to {kg_unprocessed_mt}")
    mt.write(kg_unprocessed_mt, overwrite=True)


# TODO: This function as a full copy form the step 2/2-prune-related_samples. Need to extract it to a separate module
def run_pc_relate(pruned_mt_file: str, relatedness_ht_file: str, samples_to_remove_file: str, scores_file: str) -> None:
    """
    Runs PC relate on pruned MT
    :param str pruned_mt_file: ld pruned MT file
    :param str relatedness_ht_file: relatedness ht file
    :param str samples_to_remove_file: file samples to remove ht is written to
    :param str scores_file: file scores ht is written to
    """
    print("=== Running PC relate")
    pruned_mt = hl.read_matrix_table(pruned_mt_file)
    eig, scores, _ = hl.hwe_normalized_pca(pruned_mt.GT, k=10, compute_loadings=False)
    scores.write(scores_file, overwrite=True)

    print("=== Calculating relatedness (this step usually takes a while)")
    relatedness_ht = hl.pc_relate(
        pruned_mt.GT,
        min_individual_maf=0.05,
        scores_expr=scores[pruned_mt.col_key].scores,
        block_size=4096,
        min_kinship=0.05,
        statistics="kin2",
    )
    relatedness_ht.write(relatedness_ht_file, overwrite=True)
    # prune individuals to be left with unrelated - creates a table containing one column - samples to remove
    pairs = relatedness_ht.filter(relatedness_ht["kin"] > 0.125)
    related_samples_to_remove = hl.maximal_independent_set(pairs.i, pairs.j, keep=False)
    related_samples_to_remove.write(samples_to_remove_file, overwrite=True)


def get_options() -> Any:
    """
    Get options from the command line
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--kg-to-mt", help="Convert 1kg data to matrixtable", action="store_true")
    parser.add_argument("--kg-filter-and-prune", help="Prune related variants from 1kg matrix", action="store_true")
    parser.add_argument("--kg-pc-relate", help="Run PC relate for 1KG ", action="store_true")
    parser.add_argument("--kg-remove-related-samples", help="Run PC relate for 1KG", action="store_true")
    parser.add_argument("-a", "--all", help="Run All steps", action="store_true")
    args = parser.parse_args()
    return args


def main() -> None:
    args = get_options()
    # set up input variables
    config = parse_config()
    tmp_dir = config["general"]["tmp_dir"]

    # initialise hail
    _ = hail_utils.init_hl(tmp_dir)

    step_conf = config["step1"]["create_1kg_mt"]
    vcf_indir = path_spark(step_conf["indir"])
    kg_unprocessed_mt = path_spark(step_conf["kg_unprocessed"])
    if args.kg_to_mt or args.all:
        create_1kg_mt(vcf_indir, kg_unprocessed_mt)

    pruned_kg_file = path_spark(step_conf["pruned_kg_file"])
    long_range_ld_file = path_spark(step_conf["long_range_ld_file"])
    call_rate_threshold = float(step_conf["call_rate_threshold"])
    af_threshold = float(step_conf["af_threshold"])
    hwe_threshold = float(step_conf["hwe_threshold"])

    if args.kg_filter_and_prune or args.all:
        # prunning of the linked Variants
        kg_mt = hl.read_matrix_table(kg_unprocessed_mt)
        kg_mt_filtered = filtering.filter_matrix_for_ldprune(
            kg_mt, long_range_ld_file, call_rate_threshold, af_threshold, hwe_threshold
        )

        # TODO: this par from the the file 2/2-3a-merge-and-prune. Need to extract to a separate function
        # prunning some part of the chromosome,  the linked Variantion regions that are related to each other
        pruned_kg_ht = hl.ld_prune(kg_mt_filtered.GT, r2=0.2)
        pruned_kg_mt = kg_mt_filtered.filter_rows(hl.is_defined(pruned_kg_ht[kg_mt_filtered.row_key]))
        pruned_kg_mt = pruned_kg_mt.select_entries(
            GT=hl.unphased_diploid_gt_index_call(pruned_kg_mt.GT.n_alt_alleles())
        )
        # saving matrix
        pruned_kg_mt.write(pruned_kg_file, overwrite=True)

    # TODO: this part is from the file 2-prune-related-samples. NExx to extract to separate function
    relatedness_ht_file = path_spark(step_conf["relatedness_ht_file"])
    samples_to_remove_file = path_spark(step_conf["samples_to_remove_file"])
    scores_file = path_spark(step_conf["scores_file"])
    if args.kg_pc_relate or args.all:
        run_pc_relate(pruned_kg_file, relatedness_ht_file, samples_to_remove_file, scores_file)

    kg_mt_file = path_spark(step_conf["kg_out_mt"])
    if args.kg_remove_related_samples or args.all:
        # TODO: This part of code is from the script 2/2 fucntion run_population_pca().
        kg_mt = hl.read_matrix_table(kg_unprocessed_mt)
        variants, samples = kg_mt.count()
        print(f"=== Loaded form initial table: {samples} samples, {variants} variants.")
        print("=== Removing related samples")
        related_samples_to_remove = hl.read_table(samples_to_remove_file)
        kg_mt_remove_related = kg_mt.filter_cols(hl.is_defined(related_samples_to_remove[kg_mt.col_key]), keep=False)
        variants, samples = kg_mt_remove_related.count()
        print(f"=== Remains after removing related samples: {samples} samples, {variants} variants.")
        kg_mt_remove_related.write(kg_mt_file, overwrite=True)

    # hail_utils.stop_hl(sc)


if __name__ == "__main__":
    main()
