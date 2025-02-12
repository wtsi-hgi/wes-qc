# compare different combinations of hard filters
import argparse
import os.path
from pathlib import Path

import hail as hl
import json
import logging
from typing import Optional, Any, Union
from utils.utils import parse_config, path_spark
from utils.utils import select_founders, collect_pedigree_samples
from wes_qc import hail_utils
import pandas as pd
import bokeh.plotting
import bokeh.layouts
import bokeh.io
import bokeh.transform
import bokeh.palettes
import bokeh.models
import numpy as np
import time
import datetime

snv_label = "snv"
indel_label = "indel"


def clean_mt(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Reduces matrixtable by cleaning all fields not required for hardfilter evaluation
    """
    mt = mt.select_entries(mt.GT, mt.HetAB, mt.DP, mt.GQ)
    mt = mt.drop(mt.assigned_pop, *mt.row_value)
    mt = mt.annotate_rows(info=hl.Struct())
    mt = mt.annotate_rows(
        type=hl.case()
        .when(hl.is_snp(mt.alleles[0], mt.alleles[1]), snv_label)
        .when(hl.is_indel(mt.alleles[0], mt.alleles[1]), indel_label)
        .default("other")
    )
    return mt


def annotate_with_rf(mt: hl.MatrixTable, rf_ht: hl.Table) -> hl.MatrixTable:
    """
    Annotate MatrixTable with TP, FP, rf_bin and rf_score
    :param hl.MatrixTable mt: Input MatrixTable
    :param str rf_htfile: Random forest ht file
    :return: hl.MatrixTable
    """
    # keep only vars with rank_id = rank
    rf_ht = rf_ht.filter(rf_ht.rank_id == "rank")

    # annotate mt with score and bin
    mt = mt.annotate_rows(info=mt.info.annotate(rf_score=rf_ht[mt.row_key].score))
    mt = mt.annotate_rows(info=mt.info.annotate(rf_bin=rf_ht[mt.row_key].bin))

    mt = mt.annotate_rows(TP=rf_ht[mt.row_key].tp)
    mt = mt.annotate_rows(FP=rf_ht[mt.row_key].fp)

    return mt


def annotate_cq(mt: hl.MatrixTable, cqfile: str) -> hl.MatrixTable:
    """
    Annotate MatrixTable with variant consequences taken from VEP annotation
    :param hl.MatrixTable mt: Input MatrixTable
    :param str cqfile: Most severe consequence annotation from VEP
    :return: hl.MatrixTable
    """
    ht = hl.import_table(path_spark(cqfile), types={"f1": "int32"}, no_header=True)
    ht = ht.rename({"f0": "chr", "f1": "pos", "f2": "rs", "f3": "ref", "f4": "alt", "f5": "consequence"})
    ht = ht.key_by(locus=hl.locus(ht.chr, ht.pos), alleles=[ht.ref, ht.alt])
    ht = ht.drop(ht.chr, ht.pos, ht.ref, ht.alt)
    mt = mt.annotate_rows(consequence=ht[mt.row_key].consequence)
    return mt


def prepare_giab_ht(giab_vcf: str, giab_cqfile: str) -> hl.Table:
    """
    Make GIAB table from the VCF file
    :param str giab_vcf: input VCF file
    :param str giab_cqfile: Vep annotation
    :return: hl.Table
    """
    logging.info("Preparing GiaB HailTable")
    mt = hl.import_vcf(path_spark(giab_vcf), force_bgz=True, reference_genome="GRCh38")
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(mt.variant_qc.n_non_ref > 0)

    ht = hl.import_table(path_spark(giab_cqfile), types={"f1": "int32"}, no_header=True)
    ht = ht.rename(
        {
            "f0": "chr",
            "f1": "pos",
            "f2": "rs",
            "f3": "ref",
            "f4": "alt",
            "f5": "consequence",
            "f6": "impacte",
            "f7": "misc",
        }
    )
    ht = ht.key_by(locus=hl.locus(ht.chr, ht.pos), alleles=[ht.ref, ht.alt])
    ht = ht.drop(ht.chr, ht.pos, ht.ref, ht.alt)

    mt = mt.annotate_rows(consequence=ht[mt.row_key].consequence)
    giab_vars = mt.rows()

    return giab_vars


def str_timedelta(delta: datetime.timedelta) -> str:
    """
    Convert a timedelta object to a human-readable string
    """
    seconds_in_hr = 60 * 60
    seconds_in_day = seconds_in_hr * 24
    seconds = int(delta.total_seconds())
    days = seconds // seconds_in_day
    hours = (seconds % seconds_in_day) / seconds_in_hr
    return f"{days} days and {hours:4.2f} hr"


def cache_filter_results(func):
    """
    Decorator to handle caching of filter combination results to JSON files.
    If a cached result exists, it will be loaded instead of recomputing.
    If not, the function will be executed and its results cached.
    """

    def wrapper(filter_name: str, json_dump_folder: str, var_type: str, **kwargs) -> dict:
        json_dump_file = os.path.join(json_dump_folder, f"{var_type}_hardfilters_{filter_name}.json")

        if os.path.exists(json_dump_file):
            with open(json_dump_file, "r") as f:
                checkpoint = json.load(f)
            print(f"--- Checkpoint data loaded from file {json_dump_file}")
            return checkpoint[filter_name]

        var_counts = func(var_type=var_type, **kwargs)

        with open(json_dump_file, "w") as f:
            json.dump({filter_name: var_counts}, f)

        print(f"--- Temporary results dumped to JSON: {json_dump_file}")
        return var_counts

    return wrapper


@cache_filter_results
def process_filter_combination(
    mt_bin: hl.MatrixTable,
    dp: int,
    gq: int,
    ab: float,
    call_rate: float,
    ht_giab_control: hl.Table,
    pedigree: hl.Pedigree,
    var_type: str,
    mtdir: str,
    giab_sample_id: Optional[str] = None,
) -> dict:
    """
    Process a single filter combination and return the variant counts and metrics.

    :param hl.MatrixTable mt_bin: Input MatrixTable filtered by bin
    :param int dp: Depth threshold
    :param int gq: Genotype quality threshold
    :param float ab: Allele balance threshold
    :param float call_rate: Call rate threshold
    :param hl.Table ht_giab_control: GIAB variants table from the external sample.
    Used as control
    :param hl.Pedigree pedigree: Hail pedigree object
    :param str var_type: Variant type (snv/indel)
    :param str mtdir: MatrixTable directory
    :param str giab_sample_id: The ID of GIAB sample in the dataset.
    If None, then no prec/recall is calculated.
    :return: Dictionary containing variant counts and metrics
    """
    # Applying the combination of hard filters
    mt_hard_path = os.path.join(mtdir, f"tmp.hard_filters_combs.{var_type}-hard.mt")
    mt_hard = apply_hard_filters(mt_bin, dp=dp, gq=gq, ab=ab, call_rate=call_rate, checkpoint_path=mt_hard_path)

    var_counts = {}
    # Counting TP and FP
    var_counts["TP"], var_counts["FP"] = count_tp_fp(mt_hard)

    # Counting transmitted/untransmitted ratio - for SNPs only
    if var_type == "snv":
        # Extracting synonymous for transmitted/unstransmitted calculation
        mt_syn = mt_hard.filter_rows(mt_hard.consequence == "synonymous_variant")
        ratio = count_trans_untrans(mt_syn, pedigree)
        var_counts["t_u_ratio"] = ratio

    # Calculating precision/recall - possible only if we have GIAB sample
    if giab_sample_id is not None:
        # Extracting GIAB sample for precision/recall calculation
        tmp_ht_giab_dataset_path = os.path.join(mtdir, "pr.ht")
        ht_giab_dataset = get_giab_sample_from_dataset(
            giab_sample_id, mt_hard, checkpoint_path=tmp_ht_giab_dataset_path
        )

        if var_type == "snv":
            prec, recall = count_precision_recall(ht_giab_control, ht_giab_dataset)
            var_counts["prec"] = prec
            var_counts["recall"] = recall
        elif var_type == "indel":
            prec, recall, prec_frameshift, recall_frameshift, prec_inframe, recall_inframe = count_prec_recall_indel(
                ht_giab_control=ht_giab_control, ht_giab_dataset=ht_giab_dataset, mtdir=mtdir
            )
            var_counts["prec"] = prec
            var_counts["recall"] = recall
            var_counts["prec_inframe"] = prec_inframe
            var_counts["recall_inframe"] = recall_inframe
            var_counts["prec_frameshift"] = prec_frameshift
            var_counts["recall_frameshift"] = recall_frameshift
    else:
        var_counts["prec"] = -1
        var_counts["recall"] = -1
        if var_type == "indel":
            var_counts["prec_frameshift"] = -1
            var_counts["prec_inframe"] = -1
            var_counts["recall_inframe"] = -1
            var_counts["prec_frameshift"] = -1
            var_counts["recall_frameshift"] = -1

    return var_counts


def get_giab_sample_from_dataset(giab_sample_id: str, mt_hard: hl.MatrixTable, checkpoint_path):
    """
    Extracts the GIAB sample table from the dataset matrixtable
    """
    mt_giab_sample = mt_hard.filter_cols(mt_hard.s == giab_sample_id)  # GIAB sample for precision/recall
    mt_giab_sample = mt_giab_sample.filter_rows(mt_giab_sample.locus.in_autosome())
    mt_giab_sample = hl.variant_qc(mt_giab_sample)
    mt_giab_sample = mt_giab_sample.filter_rows(mt_giab_sample.variant_qc.n_non_ref > 0)
    ht_giab_dataset = mt_giab_sample.rows()
    ht_giab_dataset = ht_giab_dataset.checkpoint(checkpoint_path, overwrite=True)
    return ht_giab_dataset


def filter_and_count(
    mt: hl.MatrixTable,
    ht_giab_control: hl.Table,
    pedigree: hl.Pedigree,
    mtdir: str,
    var_type: str,
    hardfilter_combinations: dict[str, Union[int, float]],
    giab_sample_id: Optional[str],
    json_dump_folder: str,
    **kwargs,
) -> dict:
    """
    Filter MT by various bins followed by genotype GQ and calculate % of FP and TP remaining for each bin
    """

    snv_bins: int = hardfilter_combinations["snp_bins"]
    indel_bins: int = hardfilter_combinations["indel_bins"]
    gq_vals: int = hardfilter_combinations["gq_vals"]
    dp_vals: int = hardfilter_combinations["dp_vals"]
    ab_vals: float = hardfilter_combinations["ab_vals"]
    missing_vals: float = hardfilter_combinations["missing_vals"]

    Path(json_dump_folder).mkdir(parents=True, exist_ok=True)

    # Subsetting GIAB table to SNVs or indels, depending on the input parameters
    giab_ht_type_name = os.path.join(mtdir, f"giab_table.{var_type}.ht")
    if var_type == "snv":
        bins = snv_bins
        ht_giab_control = ht_giab_control.filter(
            hl.is_snp(ht_giab_control.alleles[0], ht_giab_control.alleles[1])
        ).checkpoint(giab_ht_type_name, overwrite=True)
    elif var_type == "indel":
        bins = indel_bins
        ht_giab_control = ht_giab_control.filter(
            hl.is_indel(ht_giab_control.alleles[0], ht_giab_control.alleles[1])
        ).checkpoint(giab_ht_type_name, overwrite=True)
    else:
        print(f"Unknown variant type: {var_type}")
        raise ValueError("The filter_and_count_ function can be called only for snv or for indels")

    start_time = datetime.datetime.now()

    results = {var_type: {}}

    sample_ids = set(mt.s.collect())
    if giab_sample_id is not None and giab_sample_id not in sample_ids:
        raise ValueError(f"=== Control sample {giab_sample_id} not found in matrixtable")

    n_step = 0
    total_steps = len(bins) * len(dp_vals) * len(gq_vals) * len(ab_vals) * len(missing_vals)
    print(f"=== Starting parameter combination evaluation: {total_steps} combinations ===")

    print(f"=== Starting evaluation for {var_type} ===")
    mt = mt.filter_rows(mt.type == var_type)
    n_steps_run = 0

    mt_bin_path_previous = os.path.join(mtdir, f"tmp.hard_filters_combs_{var_type}.bin.1.mt")
    mt_bin_path_current = os.path.join(mtdir, f"tmp.hard_filters_combs_{var_type}.bin.2.mt")

    # Making the first bin-filtered matrixtable
    # For the zero iteration, this is full initial matrixtable
    mtbin = mt
    mtbin.checkpoint(mt_bin_path_previous, overwrite=True)

    for bin in sorted(bins, reverse=True):
        bin_str = f"bin_{bin}"
        print(f"=== Processing {var_type} bin: {bin} ===")
        # Making next bin-filtered matrixtable form the previous matrixtable,
        # It works because we cycle from the most relaxed bin to the most stringent
        mt_bin = mt.filter_rows(mtbin.info.rf_bin <= bin)
        mt_bin = mt_bin.checkpoint(mt_bin_path_current, overwrite=True)
        # Swapping paths
        # The current path becomes previous, and the old previous becomes current for the next bin
        mt_bin_path_previous, mt_bin_path_current = mt_bin_path_current, mt_bin_path_previous

        for dp in dp_vals:
            dp_str = f"DP_{dp}"
            for gq in gq_vals:
                gq_str = f"GQ_{gq}"
                for ab in ab_vals:
                    ab_str = f"AB_{ab}"
                    for call_rate in missing_vals:
                        missing_str = f"missing_{call_rate}"

                        print(
                            f"--- Testing {var_type} hard filter combination: bin={bin} DP={dp} GQ={gq} AB={ab} call_rate={call_rate}"
                        )
                        logging.info(f"{dp_str} {gq_str} {ab_str} {missing_str}")
                        filter_name = "_".join([bin_str, dp_str, gq_str, ab_str, missing_str])
                        # Running evaluation for the hardfilter combination
                        var_counts = process_filter_combination(
                            mt_bin=mt_bin,
                            dp=dp,
                            gq=gq,
                            ab=ab,
                            call_rate=call_rate,
                            ht_giab_control=ht_giab_control,
                            pedigree=pedigree,
                            var_type=var_type,
                            mtdir=mtdir,
                            giab_sample_id=giab_sample_id,
                            filter_name=filter_name,
                            json_dump_folder=json_dump_folder,
                        )

                        results[var_type][filter_name] = var_counts
                        # Cleaning up Hail temporary folder to avoid exceeding inodes limit.
                        hail_utils.clear_temp_folder(hl.tmp_dir())

                        stop_time = datetime.datetime.now()
                        n_steps_run += 1

                        n_step += 1
                        print(f"--- Step {n_step}/{total_steps} completed.")
                        if n_steps_run > 0:
                            total_time = stop_time - start_time
                            print(
                                f"--- Spend {total_time.days} days ({total_time.seconds / 3600:4.3f} hr) to process {n_steps_run} steps"
                            )
                            step_avg_time = total_time / n_steps_run
                            est_time = step_avg_time * (total_steps - n_step)
                            print(f"--- Estimated to complete: {str_timedelta(est_time)}.")

    print(f"=== Calculating total TP and FP for {var_type} ===")

    (results[f"{var_type}_total_tp"], results[f"{var_type}_total_fp"]) = count_tp_fp(mt)

    return results


def count_tp_fp(mt: hl.MatrixTable) -> tuple[int, int]:
    """
    Count true positives and false positives from a Hail MatrixTable.
    Parameters:
        mt (hl.MatrixTable): Input Hail MatrixTable that contains the 'TP' and 'FP' fields
                             used for counting.
    Returns:
        tuple[int, int]: A tuple (TP, FP).
    """
    ht = mt.rows()
    value_counts = ht.aggregate(hl.struct(tp=hl.agg.count_where(ht.TP), fp=hl.agg.count_where(ht.FP)))
    return value_counts.tp, value_counts.fp


def apply_hard_filters(
    mt: hl.MatrixTable, dp: int, gq: int, ab: float, call_rate: float, checkpoint_path: str
) -> hl.MatrixTable:
    """
    Applies hard filters to a Hail MatrixTable based on genotype, depth, quality, allele balance,
    and variant call rate metrics.
    mt : hl.MatrixTable. Input Hail MatrixTable to be filtered.
    dp : int. Minimum depth threshold for filtering.
    gq : int. Minimum genotype quality threshold for filtering.
    ab : float Minimum allele balance threshold for heterozygous calls.
    call_rate : float. Minimum call rate threshold for filtering variants.

    Returns
    hl.MatrixTable. Filtered Hail MatrixTable after applying hard filters and call rate filtering.
    """
    # Filtering out entries by the hardfilter combination
    filter_condition = (mt.GT.is_het() & (mt.HetAB < ab)) | (mt.DP < dp) | (mt.GQ < gq)
    mt_tmp = mt.annotate_entries(hard_filters=hl.if_else(filter_condition, "Fail", "Pass"))
    mt_tmp = mt_tmp.filter_entries(mt_tmp.hard_filters == "Pass")

    # Filter variants by call rate
    mt_tmp = mt_tmp.annotate_rows(pass_count=hl.agg.count_where(mt_tmp.hard_filters == "Pass"))
    mt_tmp = mt_tmp.filter_rows(mt_tmp.pass_count / mt_tmp.count_cols() > call_rate)

    # remove reference-only rows
    mt_tmp = hl.variant_qc(mt_tmp)
    mt_tmp = mt_tmp.filter_rows(mt_tmp.variant_qc.n_non_ref == 0, keep=False)
    mt_tmp = mt_tmp.checkpoint(checkpoint_path, overwrite=True)

    return mt_tmp


def count_prec_recall_indel(ht_giab_control: hl.Table, ht_giab_dataset: hl.Table, mtdir: str) -> tuple:
    """
    Get precison/recall vs GIAB for indels
    In fact, calls calculate_precision_recall for different subset of indels
    :param hl.Table ht_giab_dataset: hail Table of variants in ALSPAC GIAB sample
    :param hl.Table ht_giab_control: hail Table GIAB variants
    :param str var_type: Variant type
    :param str mtdir: matrixtable directory
    :return: tuple
    """
    tmp_ht_giab_dataset_indels = os.path.join(mtdir, "tmp_pr_indel.ht")
    ht_giab_dataset_indels = ht_giab_dataset.checkpoint(tmp_ht_giab_dataset_indels, overwrite=True)

    # Regular precision-recall for the full set of indels
    p, r = count_precision_recall(ht_giab_control, ht_giab_dataset_indels)

    # Precision/recall for FrameShift indels
    ht_giab_control_frameshift = ht_giab_control.filter(ht_giab_control.consequence == "frameshift_variant")
    ht_giab_dataset_frameshift = ht_giab_dataset_indels.filter(
        ht_giab_dataset_indels.consequence == "frameshift_variant"
    )
    p_f, r_f = count_precision_recall(ht_giab_control_frameshift, ht_giab_dataset_frameshift)

    # Precision/recall for In-Frame indels
    inframe_cqs = ["inframe_deletion", "inframe_insertion"]
    ht_giab_conrol_inframe = ht_giab_control.filter(hl.literal(inframe_cqs).contains(ht_giab_control.consequence))
    ht_giab_dataset_inframe = ht_giab_dataset_indels.filter(
        hl.literal(inframe_cqs).contains(ht_giab_dataset_indels.consequence)
    )
    p_if, r_if = count_precision_recall(ht_giab_conrol_inframe, ht_giab_dataset_inframe)

    return p, r, p_f, r_f, p_if, r_if


def count_precision_recall(ht_control: hl.Table, ht_dataset: hl.Table) -> tuple[float, float]:
    """
    Calculates precision/recall for the gived control and test Hail tables.
    :param hl.Table ht_control: Control set
    :paran hl.Table ht_test: Test set
    :return tuple:
    """
    vars_in_both = ht_control.semi_join(ht_dataset)
    control_only = ht_control.anti_join(ht_dataset)
    test_only = ht_dataset.anti_join(ht_control)
    tp = vars_in_both.count()
    fn = control_only.count()
    fp = test_only.count()

    precision = tp / (tp + fp) if tp + fp > 0 else 0.0
    recall = tp / (tp + fn) if tp + fn > 0 else 0.0

    return precision, recall


def count_trans_untrans(mt_syn: hl.MatrixTable, pedigree: hl.Pedigree) -> float:
    """
    Calculates transmitted/untransmitted ratio for the synonymous matrixtable
    :param hl.MatrixTable mt_syn: matrixtable, containing only synonymous variants
    :param hl.Pedigree pedigree: Hail Pedigree
    :param str mtdir: matrixtable directory
    :return float:
    """
    # list of samples in trios
    sample_list = collect_pedigree_samples(pedigree)

    # restrict to samples in trios, annotate with AC and filter to AC == 1 in parents
    mt2 = mt_syn.filter_cols(hl.set(sample_list).contains(mt_syn.s))

    founders = select_founders(pedigree)
    mt_founders = mt2.filter_cols(hl.set(founders).contains(mt2.s))
    mt_founders = hl.variant_qc(mt_founders, name="varqc_founders")

    mt2 = mt2.annotate_rows(varqc_trios=hl.Struct(AC=mt_founders.index_rows(mt2.row_key).varqc_founders.AC))

    # split to potentially transmitted/untransmitted
    trans_mt = mt2.filter_rows(mt2.varqc_trios.AC[1] == 1)

    # run tdt function for potential trans and untrans
    tdt_ht = hl.transmission_disequilibrium_test(trans_mt, pedigree)

    # Collecting stats in a single pass
    t_u_stats = tdt_ht.aggregate(hl.struct(trans=hl.agg.sum(tdt_ht.t), untrans=hl.agg.sum(tdt_ht.u)))

    if t_u_stats.untrans > 0:
        ratio = t_u_stats.trans / t_u_stats.untrans
    else:
        ratio = -1

    return ratio


def write_snv_filter_metrics(results: dict, outfile: str):
    """
    Write SNV filtering metrics to a tab-delimited file, including true/false positive rates,
    precision/recall, and transmission/untransmission ratio.

    :param dict results: Dictionary containing filtering results and metrics for each SNV filter combination
    :param str outfile: Path to output TSV file
    """
    header = ["filter", "bin", "DP", "GQ", "AB", "call_rate", "TP", "FP", "t_u_ratio", "precision", "recall"]

    with open(outfile, "w") as o:
        o.write("\t".join(header))
        o.write("\n")

        for var_f in results["snv"].keys():
            bin_val, dp, gq, ab, call_rate = parse_hard_filter_values(var_f)
            tp = str((results["snv"][var_f]["TP"] / results["snv_total_tp"]) * 100)
            fp = str((results["snv"][var_f]["FP"] / results["snv_total_fp"]) * 100)
            tu = str(results["snv"][var_f]["t_u_ratio"])
            p = str(results["snv"][var_f].get("prec", ""))
            r = str(results["snv"][var_f].get("recall", ""))

            outline = [var_f, bin_val, dp, gq, ab, call_rate, tp, fp, tu, p, r]
            o.write("\t".join(outline))
            o.write("\n")


def write_indel_filter_metrics(results: dict, outfile: str):
    """
    Write indel filtering metrics to a tab-delimited file, including true/false positive rates,
    precision/recall for all indels, frameshift indels, and inframe indels.

    :param dict results: Dictionary containing filtering results and metrics for each indel filter combination
    :param str outfile: Path to output TSV file
    """
    header = [
        "filter",
        "bin",
        "DP",
        "GQ",
        "AB",
        "call_rate",
        "TP",
        "FP",
        "precision",
        "recall",
        "precision_frameshift",
        "recall_frameshift",
        "precision_inframe",
        "recall_inframe",
    ]

    with open(outfile, "w") as o:
        o.write("\t".join(header))
        o.write("\n")

        for var_f in results["indel"].keys():
            bin_val, dp, gq, ab, call_rate = parse_hard_filter_values(var_f)
            tp = str((results["indel"][var_f]["TP"] / results["indel_total_tp"]) * 100)
            fp = str((results["indel"][var_f]["FP"] / results["indel_total_fp"]) * 100)
            p = str(results["indel"][var_f].get("prec", ""))
            r = str(results["indel"][var_f].get("recall", ""))
            p_f = str(results["indel"][var_f].get("prec_frameshift", ""))
            r_f = str(results["indel"][var_f].get("recall_frameshift", ""))
            p_if = str(results["indel"][var_f].get("prec_inframe", ""))
            r_if = str(results["indel"][var_f].get("recall_inframe", ""))

            outline = [var_f, bin_val, dp, gq, ab, call_rate, tp, fp, p, r, p_f, r_f, p_if, r_if]
            o.write("\t".join(outline))
            o.write("\n")


def parse_hard_filter_values(filter_string: str) -> tuple:
    """
    Get hard filter values from a filter string
    :param str filter_string: Filter string
    :return tuple: Tuple of hard filter values, they are bin, DP, GQ, AB, call_rate
    """
    values = filter_string.split("_")
    return values[1], values[3], values[5], values[7], values[9]


def plot_hard_filter_combinations(df: pd.DataFrame, x: str, y: str, outfile: str):
    """
    Create an interactive scatter plot of hard filter combinations with enhanced controls

    :param pd.DataFrame df: Input dataframe with filter metrics
    :param str x: Column name for x-axis
    :param str y: Column name for y-axis
    :param str outfile: Path to save the plot
    """
    # Create the data source
    # Convert DP to string for factor mapping
    df["DP"] = df["DP"].astype(str)
    source = bokeh.models.ColumnDataSource(df)
    filtered_source = bokeh.models.ColumnDataSource(data=source.data)
    margin_x = abs(df[x].max() - df[x].min()) * 0.05
    max_x = df[x].max() + margin_x
    min_x = df[x].min() - margin_x
    margin_y = abs(df[y].max() - df[y].min()) * 0.05
    max_y = df[y].max() + margin_y
    min_y = df[y].min() - margin_y

    # Define markers for different DP values
    MARKERS = ["circle", "diamond", "triangle", "square", "star", "plus"]
    DPs = [str(i) for i in sorted(df["DP"].unique().tolist())]

    tooltips = [
        ("--", "--"),
        ("Bin", "@bin"),
        (x, f"@{x}"),
        (y, f"@{y}"),
        ("DP", "@DP"),
        ("GQ", "@GQ"),
        ("AB", "@AB"),
        ("Call Rate", "@call_rate"),
    ]
    hover = bokeh.models.HoverTool(tooltips=tooltips)
    plot = bokeh.plotting.figure(
        title=(" ").join([y, "v", x]),
        height=800,
        width=1400,
        x_axis_label=x,
        y_axis_label=y,
        x_range=(min_x, max_x),
        y_range=(min_y, max_y),
        tools=[hover, "pan", "box_zoom", "wheel_zoom", "reset"],
    )

    bin_min = df["bin"].min()
    bin_max = df["bin"].max()
    unique_bins = sorted(df["bin"].unique())
    num_bins = len(unique_bins)
    interval = (bin_max - bin_min) / num_bins

    # Define available palettes
    palette_dict = {
        "Rainbow": bokeh.palettes.TolRainbow23[:: max(1, len(bokeh.palettes.TolRainbow23) // num_bins)][:num_bins],
        "Turbo": bokeh.palettes.Turbo256[:: max(1, len(bokeh.palettes.Turbo256) // num_bins)][:num_bins],
    }

    # Create palette selector
    palette_select = bokeh.models.Select(title="Color Palette", value="Turbo", options=list(palette_dict.keys()))

    # Initialize color mapper with default palette
    color_mapper = bokeh.models.LinearColorMapper(palette=palette_dict["Turbo"], low=bin_min, high=bin_max)

    # Add scatter points with continuous color mapping and marker mapping
    plot.scatter(
        x,
        y,
        source=filtered_source,
        size=14,
        alpha=0.6,
        color={"field": "bin", "transform": color_mapper},
        marker=bokeh.transform.factor_mark("DP", MARKERS, DPs),
        legend_group="DP",
    )

    # Configure legend
    plot.legend.click_policy = "hide"
    plot.legend.location = "top_left"
    plot.legend.title = "DP"
    plot.legend.background_fill_alpha = 0.7

    # Add color bar with custom ticks
    ticks = np.linspace(bin_min + interval / 2, bin_max - interval / 2, num_bins)
    color_bar = bokeh.models.ColorBar(
        color_mapper=color_mapper,
        label_standoff=12,
        title="Bin",
        location=(0, 0),
        orientation="vertical",
        ticker=bokeh.models.FixedTicker(ticks=ticks),
        formatter=bokeh.models.CustomJSTickFormatter(
            code="""
                            return bins[index].toString();
                        """,
            args={"bins": list([int(i) for i in unique_bins])},
        ),
    )

    plot.add_layout(color_bar, "right")

    # Create checkboxes for filtering
    unique_dp = [str(i) for i in sorted(df["DP"].unique().tolist())]
    unique_gq = [str(i) for i in sorted(df["GQ"].unique().tolist())]
    unique_ab = [f"{i:.2f}" for i in sorted(df["AB"].unique().tolist())]
    unique_cr = [f"{i:.2f}" for i in sorted(df["call_rate"].unique().tolist())]

    checkbox_dp = bokeh.models.CheckboxGroup(labels=unique_dp, active=list(range(len(unique_dp))), name="DP Filter")
    checkbox_gq = bokeh.models.CheckboxGroup(labels=unique_gq, active=list(range(len(unique_gq))), name="GQ Filter")
    checkbox_ab = bokeh.models.CheckboxGroup(labels=unique_ab, active=list(range(len(unique_ab))), name="AB Filter")
    checkbox_cr = bokeh.models.CheckboxGroup(
        labels=unique_cr, active=list(range(len(unique_cr))), name="Call Rate Filter"
    )

    # Create bin filter sliders
    min_bin_slider = bokeh.models.Slider(
        start=bin_min, end=bin_max, value=bin_min, step=1, title="Minimum Bin", width=200
    )
    max_bin_slider = bokeh.models.Slider(
        start=bin_min, end=bin_max, value=bin_max, step=1, title="Maximum Bin", width=200
    )

    # JavaScript callback for filtering based on all controls
    callback = bokeh.models.CustomJS(
        args=dict(
            source=source,
            filtered_source=filtered_source,
            checkbox_dp=checkbox_dp,
            checkbox_gq=checkbox_gq,
            checkbox_ab=checkbox_ab,
            checkbox_cr=checkbox_cr,
            color_mapper=color_mapper,
            palette_select=palette_select,
            palette_dict=palette_dict,
            min_bin_slider=min_bin_slider,
            max_bin_slider=max_bin_slider,
        ),
        code="""
        const active_labels_dp = checkbox_dp.active.map(i => checkbox_dp.labels[i]);
        const active_labels_gq = checkbox_gq.active.map(i => checkbox_gq.labels[i]);
        const active_labels_ab = checkbox_ab.active.map(i => checkbox_ab.labels[i]);
        const active_labels_cr = checkbox_cr.active.map(i => checkbox_cr.labels[i]);

        const data = source.data;
        const filtered_data = {};

        Object.keys(data).forEach(key => {
            filtered_data[key] = [];
        });

        for (let i = 0; i < data['DP'].length; i++) {
            const match_gq = active_labels_gq.includes(String(data['GQ'][i]));
            const match_ab = active_labels_ab.includes(data['AB'][i].toFixed(2));
            const match_cr = active_labels_cr.includes(data['call_rate'][i].toFixed(2));
            const match_dp = active_labels_dp.includes(String(data['DP'][i]));
            const match_bin = data['bin'][i] >= min_bin_slider.value &&
                            data['bin'][i] <= max_bin_slider.value;

            if (match_dp && match_gq && match_ab && match_cr && match_bin) {
                Object.keys(data).forEach(key => {
                    filtered_data[key].push(data[key][i]);
                });
            }
        }

        filtered_source.data = filtered_data;
        filtered_source.change.emit();
        """,
    )

    # Add palette change callback
    palette_callback = bokeh.models.CustomJS(
        args=dict(color_mapper=color_mapper, palette_dict=palette_dict),
        code="""
        const new_palette = palette_dict[cb_obj.value];
        color_mapper.palette = new_palette;
        """,
    )

    # Connect all callbacks
    min_bin_slider.js_on_change("value", callback)
    max_bin_slider.js_on_change("value", callback)
    checkbox_dp.js_on_change("active", callback)
    checkbox_gq.js_on_change("active", callback)
    checkbox_ab.js_on_change("active", callback)
    checkbox_cr.js_on_change("active", callback)
    palette_select.js_on_change("value", palette_callback)

    # Create layout with all controls
    controls = bokeh.layouts.column(
        bokeh.layouts.row(
            bokeh.layouts.column(bokeh.models.Div(text="<b>Genotype Depth (DP)</b>"), checkbox_dp, width=100),
            bokeh.layouts.column(bokeh.models.Div(text="<b>Genotype Quality (GQ)</b>"), checkbox_gq, width=100),
        ),
        bokeh.layouts.row(
            bokeh.layouts.column(bokeh.models.Div(text="<b>Allele Balance (AB)</b>"), checkbox_ab, width=100),
            bokeh.layouts.column(bokeh.models.Div(text="<b>Call Rate</b>"), checkbox_cr, width=100),
        ),
        bokeh.models.Div(text="<b>Bin Range Filter</b>"),
        bokeh.layouts.row(bokeh.layouts.column(min_bin_slider, max_bin_slider, width=200)),
        bokeh.layouts.row(bokeh.layouts.column(palette_select, width=200)),
    )
    layout = bokeh.layouts.row(controls, plot)

    # Save to file
    bokeh.io.output_file(outfile)
    bokeh.io.save(layout)


def get_options() -> Any:
    """
    Get options from the command line
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepare", help="Prepare all required matrixtables", action="store_true")
    parser.add_argument("--evaluate-snv", help="Run hardfilter evaluation for SNV", action="store_true")
    parser.add_argument("--evaluate-indel", help="Run hardfilter evaluation for InDel", action="store_true")
    parser.add_argument("--plot", help="Plot hard filter combinations", action="store_true")
    parser.add_argument("--all", help="Run All steps", action="store_true")
    args = parser.parse_args()
    return args


def main():
    # = STEP SETUP = #
    config = parse_config()
    args = get_options()
    if args.all:
        args.prepare = True
        args.evaluate_snv = True
        args.evaluate_indel = True
        args.plot = True

    tmp_dir = config["general"]["tmp_dir"]

    # = STEP PARAMETERS = #
    model_id = config["general"]["rf_model_id"]
    mtdir = config["general"]["matrixtables_dir"]
    rf_dir = config["general"]["var_qc_rf_dir"]

    hardfilter_evaluate_workdir = os.path.join(mtdir, model_id)

    # = STEP DEPENDENCIES = #
    # GIAB sample to compare
    giab_vcf = config["step4"]["evaluation"]["giab_vcf"]
    giab_cqfile = config["step4"]["evaluation"]["giab_cqfile"]

    # Files from VariantQC
    mtfile = config["step3"]["split_multi_and_var_qc"]["varqc_mtoutfile_split"]
    cqfile = config["step3"]["annotate_mt_with_cq_rf_score_and_bin"]["cqfile"]
    pedfile = config["step3"]["pedfile"]
    rf_htfile = os.path.join(rf_dir, model_id, "_gnomad_score_binning_tmp.ht")

    # = STEP OUTPUTS = #
    mt_annot_path = os.path.join(hardfilter_evaluate_workdir, "tmp.hard_filters_combs.mt")
    outfile_snv = config["step4"]["evaluation"]["snp_tsv"]
    outfile_indel = config["step4"]["evaluation"]["indel_tsv"]
    giab_ht_file = config["step4"]["evaluation"]["giab_ht_file"]

    # = STEP LOGIC = #
    hail_utils.init_hl(tmp_dir)

    if args.prepare:
        print("=== Preparing control GIAB sample ===")
        ht_giab_control = prepare_giab_ht(giab_vcf, giab_cqfile)
        ht_giab_control.write(path_spark(giab_ht_file), overwrite=True)

        print("=== Annotating matrix table with RF bins ===")
        mt = hl.read_matrix_table(path_spark(mtfile))
        mt = clean_mt(mt)  # Remove all information not required for hard filter evaluation

        rf_ht = hl.read_table(path_spark(rf_htfile))
        mt_annot = annotate_with_rf(mt, rf_ht)
        mt_annot = annotate_cq(mt_annot, cqfile)
        mt_annot.write(path_spark(mt_annot_path), overwrite=True)

    if args.evaluate_snv:
        print("=== Calculating hard filter evaluation for SNV ===")
        start_time = time.time()
        mt = hl.read_matrix_table(path_spark(mt_annot_path))
        ht_giab_control = hl.read_table(path_spark(giab_ht_file)) if giab_vcf is not None else None
        pedigree = hl.Pedigree.read(path_spark(pedfile))
        results = filter_and_count(
            mt,
            ht_giab_control,
            pedigree,
            mtdir=path_spark(hardfilter_evaluate_workdir),
            var_type="snv",
            **config["step4"]["evaluation"],
        )

        os.makedirs(os.path.dirname(outfile_snv), exist_ok=True)
        write_snv_filter_metrics(results, outfile_snv)
        end_time = time.time()
        print(f"=== SNV evaluation competed successfully.\nExecution time: {end_time - start_time:.2f} seconds ===")

    if args.evaluate_indel:
        print("=== Calculating hard filter evaluation for InDels ===")
        start_time = time.time()
        mt = hl.read_matrix_table(path_spark(mt_annot_path))
        ht_giab_control = hl.read_table(path_spark(giab_ht_file))
        pedigree = hl.Pedigree.read(path_spark(pedfile))
        results = filter_and_count(
            mt,
            ht_giab_control,
            pedigree,
            mtdir=path_spark(hardfilter_evaluate_workdir),
            var_type="indel",
            **config["step4"]["evaluation"],
        )

        os.makedirs(os.path.dirname(outfile_indel), exist_ok=True)
        write_indel_filter_metrics(results, outfile_indel)
        end_time = time.time()
        print(f"=== INDEL evaluation competed successfully.\nExecution time: {end_time - start_time:.2f} seconds ===")

    if args.plot:
        print("=== Plotting hard filter combinations ===")
        df_snv = pd.read_csv(outfile_snv, sep="\t")
        df_indel = pd.read_csv(outfile_indel, sep="\t")
        # Plot hard filter combinations
        for k, v in config["step4"]["plot"].items():
            type_, x, y = k.split("-")
            df = df_snv if type_ == "snv" else df_indel
            plot_hard_filter_combinations(df, x, y, v)
        print("=== Plotting hard filter combinations completed successfully ===")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
