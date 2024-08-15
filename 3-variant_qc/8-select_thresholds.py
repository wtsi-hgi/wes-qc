# select thresholds for random forest - for a range of thresholds for SNPs and indels
# calculate how many TPs are retained and how many FPs are included
import os
import re
from typing import Any

import hail as hl
import argparse
from utils.utils import parse_config

from wes_qc import hail_utils


def get_options() -> Any:
    """
    Get options from the command line
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--snv", required=True, help="Proposed SNV threshold. A singe value or a range")
    parser.add_argument("--indel", required=True, help="Proposed Indel threshold. A singe value or a range")
    args = parser.parse_args()
    return args


def analyse_thresholds(htfile: str, snv_threshold: str, indel_threshold: str) -> None:
    """
    For suggested thresholds calculate what proportion of TP variants are included and what propeortion of FPs
    :param str htfile: ht file annotated with rank, bin, tp and fp
    :param int snv_threshold: proposed threshold for SNVs
    :param int indel_threshold: proposed threshold for indels
    """
    # for each threshold also look at 3 bins on either side (unless <1 or >101)
    snv_t = get_threshold_range(snv_threshold)
    indel_t = get_threshold_range(indel_threshold)

    ht = hl.read_table(htfile)
    # keep only those with rank_id = rank
    ht = ht.filter(ht.rank_id == "rank")
    # split by snv and indel
    snv_ht = ht.filter(ht.allele_type == "snv")
    indel_ht = ht.filter(ht.allele_type == "snv", keep=False)

    # calculate % TP kept and % FP kept at each threshold
    snv_result = get_vars_kept(snv_ht, snv_t)
    indel_result = get_vars_kept(indel_ht, indel_t)
    # print results
    print("SNV percentage of variants kept at each threshold:")
    print_results(snv_result)
    print("Indel percentage of variants kept at each threshold:")
    print_results(indel_result)


def print_results(result: dict[int, Any]) -> None:
    """
    Print percentages of variants kept at each threshold
    :param dict result: Results to print
    """
    for t in result.keys():
        print(
            "bin "
            + str(t)
            + ": "
            + "{0:.3f}".format(result[t]["tp_pc"])
            + " percent TPs remain, "
            + "{0:.3f}".format(result[t]["fp_pc"])
            + " percent FPs remain"
        )


def get_vars_kept(ht: hl.Table, t_list: list[int]) -> dict[int, Any]:
    """
    find proportions fo TP and FPs kept in ht for each threshold in t_list
    :param hl.Table ht: input Hail Table
    :param list t_list: List of thresholds
    """
    results: dict[int, Any] = {}
    ht_TP = ht.filter(ht.tp == True)
    ht_FP = ht.filter(ht.fp == True)
    tp_total = ht_TP.count()
    fp_total = ht_FP.count()
    for t in t_list:
        results[t] = {}
        ht_tp_pass = ht_TP.filter(ht_TP.bin <= t)
        tp_pass = ht_tp_pass.count()
        results[t]["tp_pc"] = 100 * (tp_pass / tp_total)
        ht_fp_pass = ht_FP.filter(ht_FP.bin <= t)
        fp_pass = ht_fp_pass.count()
        results[t]["fp_pc"] = 100 * (fp_pass / fp_total)

    return results


def get_threshold_range(threshold_str: str) -> list[int]:
    """
    Given a single threshold value convert to a range +/- 8  (unless <1 or >101)
    Given a range notation generates the range
    :param threshold: single threshold value or a list
    """

    if re.match(r"^\d+$", threshold_str):
        threshold = int(threshold_str)
        min_t = threshold - 8
        max_t = threshold + 8
    elif match := re.match(r"^(\d+)-(\d+)$", threshold_str):
        min_t = int(match.group(1))
        max_t = int(match.group(2))

    t_list = list(range(max(1, min_t), min(max_t + 1, 101)))

    return t_list


def main() -> None:
    # set up
    args = get_options()
    inputs = parse_config()

    data_root: str = inputs["data_root"]
    dataset_name: str = inputs["dataset_name"]
    rf_dir = os.path.join(data_root, inputs["var_qc_rf_dir"])
    runhash = inputs["runhash"]

    # initialise hail
    sc = hail_utils.init_hl(inputs["tmp_dir"])

    htfile = os.path.join(rf_dir, runhash, "_gnomad_score_binning_tmp.ht")
    analyse_thresholds("file://" + htfile, args.snv, args.indel)

    hail_utils.stop_hl(sc)


if __name__ == "__main__":
    main()
