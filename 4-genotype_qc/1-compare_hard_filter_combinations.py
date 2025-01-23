# compare different combinations of hard filters
import argparse
import os.path
import hail as hl
import json
import logging
from typing import Optional, Any
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

snv_label = "snv"
indel_label = "indel"


def clean_mt(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    reduces matrixtables by cleaning all non-nesessary information
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


def annotate_with_rf(mt: hl.MatrixTable, rf_htfile: str) -> hl.MatrixTable:
    """
    Annotate MatrixTable with TP, FP, rf_bin and rf_score
    :param hl.MatrixTable mt: Input MatrixTable
    :param str rf_htfile: Random forest ht file
    :return: hl.MatrixTable
    """
    rf_ht = hl.read_table(rf_htfile)

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
    Annotate MatrixTable with consequence
    :param hl.MatrixTable mt: Input MatrixTable
    :param str cqfile: Most severe consequence annotation from VEP
    :return: hl.MatrixTable
    """
    ht = hl.import_table(cqfile, types={"f1": "int32"}, no_header=True)
    ht = ht.rename({"f0": "chr", "f1": "pos", "f2": "rs", "f3": "ref", "f4": "alt", "f5": "consequence"})

    ht = ht.key_by(locus=hl.locus(ht.chr, ht.pos), alleles=[ht.ref, ht.alt])
    ht = ht.drop(ht.chr, ht.pos, ht.ref, ht.alt)

    mt = mt.annotate_rows(consequence=ht[mt.row_key].consequence)

    return mt


def prepare_giab_ht(giab_vcf: str, giab_cqfile: str, mtdir: str) -> hl.Table:
    """
    Get GIAB ht from vcf file
    :param str giab_vcf: path of input VCF file
    :param str giab_cqfile: path of vep annotation
    :param str mtdir: MatrixTable directory
    :return: hl.Table
    """
    logging.info("Preparing GiaB HailTable")
    mt = hl.import_vcf(giab_vcf, force_bgz=True, reference_genome="GRCh38")
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(mt.variant_qc.n_non_ref > 0)

    ht = hl.import_table(giab_cqfile, types={"f1": "int32"}, no_header=True)
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
    tmphtg = path_spark(os.path.join(mtdir, "tmphtgx.ht"))
    giab_vars = giab_vars.checkpoint(tmphtg, overwrite=True)

    return giab_vars


def filter_and_count(mt_path: str, ht_giab: hl.Table, pedfile: str, mtdir: str, conf: dict) -> dict:
    """
    Filter MT by various bins followed by genotype GQ and calculate % of FP and TP remaining for each bin
    :param hl.Table ht_giab: GIAB variants
    :param str mtdir: matrixtable directory
    :param mt_path: path to matrixtable
    :param pedfile: path to trio file
    :param config: config object
    :return: dict
    """
    results = {"snv": {}, "indel": {}}
    giab_sample = conf["giab_sample"]

    mt = hl.read_matrix_table(mt_path)
    pedigree = hl.Pedigree.read(pedfile)

    mt_snv = mt.filter_rows(mt.type == snv_label)
    mt_snp_path = os.path.join(mtdir, "tmp.hard_filters_combs.snp.mt")
    mt_snv = mt_snv.checkpoint(mt_snp_path, overwrite=True)

    print("=== Starting evaluation for SNVs ===")
    snv_mt_tp, snv_mt_fp, _, _ = filter_mts(mt_snv, mtdir=mtdir, giab_sample=None)
    results["snv_total_tp"] = snv_mt_tp.count_rows()
    results["snv_total_fp"] = snv_mt_fp.count_rows()

    snp_bins = conf["snp_bins"]
    indel_bins = conf["indel_bins"]
    gq_vals = conf["gq_vals"]
    dp_vals = conf["dp_vals"]
    ab_vals = conf["ab_vals"]
    missing_vals = conf["missing_vals"]

    for bin in snp_bins:
        bin_str = f"bin_{bin}"
        logging.info(f"bin {bin}")
        print(f"--- Running evaluation for SNV bin {bin}")
        mt_snp_bin = mt_snv.filter_rows(mt_snv.info.rf_bin <= bin)

        for dp in dp_vals:
            dp_str = f"DP_{dp}"
            for gq in gq_vals:
                gq_str = f"GQ_{gq}"
                for ab in ab_vals:
                    ab_str = f"AB_{ab}"
                    for call_rate in missing_vals:
                        missing_str = f"missing_{call_rate}"
                        print(
                            f"--- Testing SNV hard filter combination: bin={bin} DP={dp} GQ={gq} AB={ab} call_rate={call_rate}"
                        )
                        logging.info(f"{dp_str} {gq_str} {ab_str} {missing_str}")
                        filter_name = "_".join([bin_str, dp_str, gq_str, ab_str, missing_str])

                        mt_snp_hard = apply_hard_filters(mt_snp_bin, dp=dp, gq=gq, ab=ab, call_rate=call_rate)
                        mt_snp_hard_path = os.path.join(mtdir, "tmp.hard_filters_combs.snp-hard.mt")
                        mt_snp_hard = mt_snp_hard.checkpoint(mt_snp_hard_path, overwrite=True)

                        mt_tp, mt_fp, mt_syn, mt_prec_recall = filter_mts(
                            mt_snp_hard, mtdir=mtdir, giab_sample=giab_sample
                        )

                        snp_counts = count_tp_fp_t_u(
                            mt_tp, mt_fp, mt_syn, mt_prec_recall, ht_giab, pedigree, "snv", mtdir
                        )
                        results["snv"][filter_name] = snp_counts

                        with open(conf["snp_json"], "w") as f:
                            json.dump(results, f)

    mt_indel = mt.filter_rows(mt.type == indel_label)
    mt_indel_path = os.path.join(mtdir, "tmp.hard_filters_combs.indel.mt")
    n_partitions = max(1, mt.n_partitions() // 4)
    mt_indel = mt_indel.repartition(n_partitions).checkpoint(mt_indel_path, overwrite=True)

    indel_mt_tp, indel_mt_fp, _, _ = filter_mts(mt_indel, mtdir=mtdir, giab_sample=None)
    results["indel_total_tp"] = indel_mt_tp.count_rows()
    results["indel_total_fp"] = indel_mt_fp.count_rows()

    for bin in indel_bins:
        bin_str = f"bin_{bin}"
        print(f"bin {bin}")
        mt_indel_bin = mt_indel.filter_rows(mt_indel.info.rf_bin <= bin)

        for dp in dp_vals:
            dp_str = f"DP_{dp}"
            for gq in gq_vals:
                gq_str = f"GQ_{gq}"
                for ab in ab_vals:
                    ab_str = f"AB_{ab}"
                    for call_rate in missing_vals:
                        missing_str = f"missing_{call_rate}"
                        print(
                            f"--- Testing InDel hard filter combination: bin={bin} DP={dp} GQ={gq} AB={ab} call_rate={call_rate}"
                        )
                        filter_name = "_".join([bin_str, dp_str, gq_str, ab_str, missing_str])

                        mt_indel_hard = apply_hard_filters(mt_indel_bin, dp=dp, gq=gq, ab=ab, call_rate=call_rate)
                        mt_indel_hard_path = os.path.join(mtdir, "tmp.hard_filters_combs.indel-hard.mt")
                        mt_indel_hard = mt_indel_hard.checkpoint(mt_indel_hard_path, overwrite=True)

                        mt_tp, mt_fp, mt_syn, mt_prec_recall = filter_mts(
                            mt_indel_hard, mtdir=mtdir, giab_sample=giab_sample
                        )

                        indel_counts = count_tp_fp_t_u(
                            mt_tp, mt_fp, mt_syn, mt_prec_recall, ht_giab, pedigree, "indel", mtdir
                        )
                        results["indel"][filter_name] = indel_counts

                        with open(conf["indel_json"], "w") as f:
                            json.dump(results, f)

    return results


def apply_hard_filters(mt: hl.MatrixTable, dp: int, gq: int, ab: float, call_rate: float) -> hl.MatrixTable:
    filter_condition = (mt.GT.is_het() & (mt.HetAB < ab)) | (mt.DP < dp) | (mt.GQ < gq)
    mt_tmp = mt.annotate_entries(hard_filters=hl.if_else(filter_condition, "Fail", "Pass"))
    mt_tmp = mt_tmp.filter_entries(mt_tmp.hard_filters == "Pass")

    mt_tmp = mt_tmp.annotate_rows(pass_count=hl.agg.count_where(mt_tmp.hard_filters == "Pass"))
    mt_tmp = mt_tmp.filter_rows(mt_tmp.pass_count / mt_tmp.count_cols() > call_rate)

    # remove unused rows
    mt_tmp = hl.variant_qc(mt_tmp)
    mt_tmp = mt_tmp.filter_rows(mt_tmp.variant_qc.n_non_ref == 0, keep=False)

    return mt_tmp


def count_tp_fp_t_u(
    mt_tp: hl.MatrixTable,
    mt_fp: hl.MatrixTable,
    mt_syn: hl.MatrixTable,
    mt_prec_recall: Optional[hl.Table],
    ht_giab: hl.Table,
    pedigree: hl.Pedigree,
    var_type: str,
    mtdir: str,
):
    """
    Filter mt by each rf bin in a list, then genotype hard filters and count remaining TP and P variants
    :param hl.MatrixTable mt_tp: Input TP MatrixTable
    :param hl.MatrixTable mt_fp: Input FP MatrixTable
    :param hl.MatrixTable mt_syn: Input synonymous MatrixTable
    :param hl.MatrixTable mt_prec_recall: mt from GIAB sample for precision/recall
    :param hl.Table ht_giab: GIAB variants
    :param hl.Pedigree pedigree: hail pedigree object
    :param  str var_type: variant type (snv/indel)
    :param str mtdir: matrixtable directory
    :return: Dict containing bin and remaning TP/FP count
    """
    results = {}

    logging.debug("Executing count_tp_fp_t_u")

    counts = count_tp_fp(mt_tp, mt_fp)
    results["TP"] = counts[0]
    results["FP"] = counts[1]

    if var_type == "snv":
        ratio = get_trans_untrans(mt_syn, pedigree, mtdir)
        results["t_u_ratio"] = ratio

    if mt_prec_recall is not None:
        ht_prec_recall = mt_prec_recall.rows()

        if var_type == "snv":
            prec, recall = get_prec_recall(ht_prec_recall, ht_giab, "snv", mtdir)
            results["prec"] = prec
            results["recall"] = recall
        elif var_type == "indel":
            prec, recall, prec_frameshift, recall_frameshift, prec_inframe, recall_inframe = get_prec_recall(
                ht_prec_recall, ht_giab, "indel", mtdir
            )
            results["prec"] = prec
            results["recall"] = recall
            results["prec_inframe"] = prec_inframe
            results["recall_inframe"] = recall_inframe
            results["prec_frameshift"] = prec_frameshift
            results["recall_frameshift"] = recall_frameshift
        else:
            raise ValueError(
                f"An incorrect function call: variant type {var_type} not supported for TP/FP and prec/recall calculations."
            )
    else:
        results["prec"] = -1
        results["recall"] = -1
        if var_type == "indel":
            results["prec_frameshift"] = -1
            results["prec_inframe"] = -1
            results["recall_inframe"] = -1
            results["prec_frameshift"] = -1
            results["recall_frameshift"] = -1
    return results


def get_prec_recall(ht_prec_recall: hl.Table, ht_giab: hl.Table, var_type: str, mtdir: str) -> tuple:
    """
    Get precison/recall vs GIAB
    :param hl.Table ht_prec_recall: hail Table of variants in ALSPAC GIAB sample
    :param hl.Table ht_giab: hail Table GIAB variants
    :param str var_type: Variant type
    :param str mtdir: matrixtable directory
    :return: tuple
    """
    tmpht1 = os.path.join(mtdir, "tmppr1x.ht")
    tmpht2 = os.path.join(mtdir, "tmppr2x.ht")
    tmpht3 = os.path.join(mtdir, "tmppr3x.ht")
    tmpht4 = os.path.join(mtdir, "tmppr4x.ht")
    tmpht5 = os.path.join(mtdir, "tmppr5x.ht")
    tmpht6 = os.path.join(mtdir, "tmppr6x.ht")

    if var_type == "snv":
        giab_snvs = ht_giab.filter(hl.is_snp(ht_giab.alleles[0], ht_giab.alleles[1]))
        alspac_snvs = ht_prec_recall.filter(hl.is_snp(ht_prec_recall.alleles[0], ht_prec_recall.alleles[1]))
        giab_snvs = giab_snvs.checkpoint(tmpht1, overwrite=True)
        alspac_snvs = alspac_snvs.checkpoint(tmpht2, overwrite=True)
        p, r = calculate_precision_recall(giab_snvs, alspac_snvs)
        return p, r

    elif var_type == "indel":
        giab_indels = ht_giab.filter(hl.is_indel(ht_giab.alleles[0], ht_giab.alleles[1]))
        alspac_indels = ht_prec_recall.filter(hl.is_indel(ht_prec_recall.alleles[0], ht_prec_recall.alleles[1]))
        giab_indels = giab_indels.checkpoint(tmpht1, overwrite=True)
        alspac_indels = alspac_indels.checkpoint(tmpht2, overwrite=True)
        p, r = calculate_precision_recall(giab_indels, alspac_indels)

        giab_frameshift = giab_indels.filter(giab_indels.consequence == "frameshift_variant")
        alspac_frameshift = alspac_indels.filter(alspac_indels.consequence == "frameshift_variant")
        giab_frameshift = giab_frameshift.checkpoint(tmpht3, overwrite=True)
        alspac_frameshift = alspac_frameshift.checkpoint(tmpht4, overwrite=True)
        p_f, r_f = calculate_precision_recall(giab_frameshift, alspac_frameshift)

        inframe_cqs = ["inframe_deletion", "inframe_insertion"]
        giab_in_frame = giab_indels.filter(hl.literal(inframe_cqs).contains(giab_indels.consequence))
        alspac_in_frame = alspac_indels.filter(hl.literal(inframe_cqs).contains(alspac_indels.consequence))
        giab_in_frame = giab_in_frame.checkpoint(tmpht5, overwrite=True)
        alspac_in_frame = alspac_in_frame.checkpoint(tmpht6, overwrite=True)
        p_if, r_if = calculate_precision_recall(giab_in_frame, alspac_in_frame)

        return p, r, p_f, r_f, p_if, r_if


def calculate_precision_recall(ht_control: hl.Table, ht_test: hl.Table) -> tuple:
    """
    Calculate orecision recall
    :param hl.Table ht_control: Control set
    :paran hl.Table ht_test: Test set
    :return tuple:
    """
    print("get intersects")
    vars_in_both = ht_control.semi_join(ht_test)
    control_only = ht_control.anti_join(ht_test)
    test_only = ht_test.anti_join(ht_control)
    print("count_vars")
    tp = vars_in_both.count()
    fn = control_only.count()
    fp = test_only.count()

    precision = tp / (tp + fp)
    recall = tp / (tp + fn)

    return precision, recall


def get_trans_untrans(mt: hl.MatrixTable, pedigree: hl.Pedigree, mtdir: str) -> float:
    """
    get transmitted/untransmitted ratio
    :param hl.MatrixTable mt: matrixtable
    :param hl.Pedigree pedigree: Hail Pedigree
    :param str mtdir: matrixtable directory
    :return float:
    """
    # list of samples in trios
    sample_list = collect_pedigree_samples(pedigree)

    # filter to synonymous
    mt_syn = mt.filter_rows(mt.consequence == "synonymous_variant")

    # restrict to samples in trios, annotate with AC and filter to AC == 1 in parents
    mt2 = mt_syn.filter_cols(hl.set(sample_list).contains(mt_syn.s))

    founders = select_founders(pedigree)
    mt_founders = mt2.filter_cols(hl.set(founders).contains(mt2.s))
    mt_founders = hl.variant_qc(mt_founders, name="varqc_founders")

    mt2 = mt2.annotate_rows(varqc_trios=hl.Struct(AC=mt_founders.index_rows(mt2.row_key).varqc_founders.AC))
    tmpmt3 = os.path.join(mtdir, "tmp1x.mt")
    mt2 = mt2.checkpoint(tmpmt3, overwrite=True)

    # split to potentially transmitted/untransmitted
    trans_mt = mt2.filter_rows(mt2.varqc_trios.AC[1] == 1)
    tmpmt5 = os.path.join(mtdir, "tmp2x.mt")
    trans_mt = trans_mt.checkpoint(tmpmt5, overwrite=True)

    # run tdt function for potential trans and untrans
    tdt_ht = hl.transmission_disequilibrium_test(trans_mt, pedigree)
    trans = tdt_ht.aggregate(hl.agg.sum(tdt_ht.t))
    untrans = tdt_ht.aggregate(hl.agg.sum(tdt_ht.u))
    if untrans > 0:
        ratio = trans / untrans
    else:
        ratio = -1

    return ratio


def count_tp_fp(mt_tp: hl.MatrixTable, mt_fp: hl.MatrixTable) -> tuple:
    """
    Count total TPs and FPs in a pair of MatrixTables
    :param hl.MatrixTable mt_tp: Input TP mt
    :param hl.MatrixTable mt_fp: Input FP mt
    :return: tuple of TP and FP counts
    """

    tp_count = mt_tp.count_rows()
    fp_count = mt_fp.count_rows()

    return tp_count, fp_count


def filter_mts(mt: hl.MatrixTable, mtdir: str, giab_sample: Optional[str] = None) -> tuple:
    """
    Split matrixtable and return tables with just TP, just FP, just synonymous
    :param hl.MatrixTable mt: Input mtfile
    :param st mtdir: matrixtable directory
    :param giab_sample: sample name of GiaB sample in cohort data
    :return: tuple of 4 hl.MatrixTable objects
    """
    mt_true = mt.filter_rows(mt.TP)  # TP variants
    mt_false = mt.filter_rows(mt.FP)  # FP variants
    mt_syn = mt.filter_rows(mt.consequence == "synonymous_variant")  # synonymous for transmitted/unstransmitted

    tmpmtt = os.path.join(mtdir, "tp.mt")
    tmpmtf = os.path.join(mtdir, "fp.mt")
    tmpmts = os.path.join(mtdir, "syn.mt")
    tmpmtpr = os.path.join(mtdir, "pr.mt")

    mt_true = mt_true.checkpoint(tmpmtt, overwrite=True)
    n_partitions = max(1, mt.n_partitions() // 10)
    mt_false = mt_false.repartition(n_partitions).checkpoint(tmpmtf, overwrite=True)
    mt_syn = mt_syn.checkpoint(tmpmts, overwrite=True)

    mt_prec_recall = None
    if giab_sample is not None:
        mt_prec_recall = mt.filter_cols(mt.s == giab_sample)  # GIAB sample for precision/recall
        mt_prec_recall = mt_prec_recall.filter_rows(mt_prec_recall.locus.in_autosome())
        mt_prec_recall = hl.variant_qc(mt_prec_recall)
        mt_prec_recall = mt_prec_recall.filter_rows(mt_prec_recall.variant_qc.n_non_ref > 0)
        mt_prec_recall = mt_prec_recall.repartition(n_partitions).checkpoint(tmpmtpr, overwrite=True)

    return mt_true, mt_false, mt_syn, mt_prec_recall


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
    parser.add_argument("--evaluate", help="Run hardfilter evaluation", action="store_true")
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
        args.evaluate = True
        args.plot = True
    tmp_dir = config["general"]["tmp_dir"]

    # = STEP PARAMETERS = #
    model_id = config["general"]["rf_model_id"]
    mtdir = config["general"]["matrixtables_dir"]
    rf_dir = config["general"]["var_qc_rf_dir"]

    wd = path_spark(os.path.join(mtdir, model_id))

    # = STEP DEPENDENCIES = #
    # GIAB sample to compare
    giab_vcf = path_spark(config["step4"]["evaluation"]["giab_vcf"])
    giab_cqfile = path_spark(config["step4"]["evaluation"]["giab_cqfile"])

    # Files from VariantQC
    mtfile = path_spark(config["step3"]["split_multi_and_var_qc"]["varqc_mtoutfile_split"])
    cqfile = path_spark(config["step3"]["annotate_mt_with_cq_rf_score_and_bin"]["cqfile"])
    pedfile = path_spark(config["step3"]["pedfile"])

    # = STEP OUTPUTS = #
    mt_annot_path = path_spark(os.path.join(wd, "tmp.hard_filters_combs.mt"))
    outfile_snv = config["step4"]["evaluation"]["snp_tsv"]
    outfile_indel = config["step4"]["evaluation"]["indel_tsv"]

    # = STEP LOGIC = #
    hail_utils.init_hl(tmp_dir)

    if args.prepare:
        print("=== Annotating matrix table with RF bins ===")
        rf_htfile = path_spark(os.path.join(rf_dir, model_id, "_gnomad_score_binning_tmp.ht"))
        mt = hl.read_matrix_table(mtfile)
        mt = clean_mt(mt)  # Remove all information not required for hard filter evaluation
        mt_annot = annotate_with_rf(mt, rf_htfile)
        mt_annot = annotate_cq(mt_annot, cqfile)
        mt_annot.write(mt_annot_path, overwrite=True)

    if args.evaluate:
        os.makedirs(os.path.dirname(outfile_snv), exist_ok=True)
        print("=== Preparing Hail table for GIAB sample ===")
        giab_ht = prepare_giab_ht(giab_vcf, giab_cqfile, mtdir)
        print("=== Calculating hard filter evaluation ===")
        results = filter_and_count(
            mt_path=mt_annot_path, ht_giab=giab_ht, pedfile=pedfile, mtdir=wd, conf=config["step4"]["evaluation"]
        )

        # Write metrics for SNVs and indels separately
        write_snv_filter_metrics(results, outfile_snv)
        write_indel_filter_metrics(results, outfile_indel)

    if args.plot:
        df_snv = pd.read_csv(outfile_snv, sep="\t")
        df_indel = pd.read_csv(outfile_indel, sep="\t")
        # Plot hard filter combinations
        for k, v in config["step4"]["plot"].items():
            type_, x, y = k.split("-")
            df = df_snv if type_ == "snv" else df_indel
            plot_hard_filter_combinations(df, x, y, v)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
