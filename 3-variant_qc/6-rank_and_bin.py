# rank and bin RF output
import hail as hl
import os.path
from typing import Optional, Dict
from pprint import pformat
from utils.utils import parse_config, path_spark
from wes_qc import hail_utils


def add_rank(
    ht: hl.Table,
    score_expr: hl.expr.NumericExpression,
    subrank_expr: Optional[Dict[str, hl.expr.BooleanExpression]] = None,
) -> hl.Table:
    """
    Adds rank based on the `score_expr`. Rank is added for snvs and indels separately.
    If one or more `subrank_expr` are provided, then subrank is added based on all sites for which the boolean expression is true.
    In addition, variant counts (snv, indel separately) is added as a global (`rank_variant_counts`).
    :param hl.Table ht: input Hail Table containing variants (with QC annotations) to be ranked
    :param hl.expr.NumericExpression score_expr: the Table annotation by which ranking should be scored
    :param subrank_expr: Any subranking to be added in the form name_of_subrank: subrank_filtering_expr
    :return: Table with rankings added
    """

    key = ht.key
    if subrank_expr is None:
        subrank_expr = {}

    temp_expr = {"_score": score_expr}
    temp_expr.update({f"_{name}": expr for name, expr in subrank_expr.items()})
    rank_ht = ht.select(**temp_expr, is_snv=hl.is_snp(ht.alleles[0], ht.alleles[1]))

    rank_ht = rank_ht.key_by("_score").persist()
    # scan_expr = {
    #     "rank": hl.cond(
    #         rank_ht.is_snv,
    #         hl.scan.count_where(rank_ht.is_snv),
    #         hl.scan.count_where(~rank_ht.is_snv),
    #     )
    # }
    scan_expr = {
        "rank": hl.if_else(
            rank_ht.is_snv,
            hl.scan.count_where(rank_ht.is_snv),
            hl.scan.count_where(~rank_ht.is_snv),
        )
    }
    # scan_expr.update(
    #     {
    #         name: hl.or_missing(
    #             rank_ht[f"_{name}"],
    #             hl.cond(
    #                 rank_ht.is_snv,
    #                 hl.scan.count_where(rank_ht.is_snv & rank_ht[f"_{name}"]),
    #                 hl.scan.count_where(~rank_ht.is_snv & rank_ht[f"_{name}"]),
    #             ),
    #         )
    #         for name in subrank_expr
    #     }
    # )
    scan_expr.update(
        {
            name: hl.or_missing(
                rank_ht[f"_{name}"],
                hl.if_else(
                    rank_ht.is_snv,
                    hl.scan.count_where(rank_ht.is_snv & rank_ht[f"_{name}"]),
                    hl.scan.count_where(~rank_ht.is_snv & rank_ht[f"_{name}"]),
                ),
            )
            for name in subrank_expr
        }
    )
    rank_ht = rank_ht.annotate(**scan_expr)

    rank_ht = rank_ht.key_by(*key).persist()
    rank_ht = rank_ht.select(*scan_expr.keys())

    ht = ht.annotate(**rank_ht[key])
    return ht


def create_binned_data_initial(
    ht: hl.Table, bin_tmp_htfile: str, truth_htfile: str, n_bins: int, config: dict
) -> hl.Table:
    """
    Create binned data from RF
    :param hl.Table ht: Input hail table
    :param str bin_tmp_htfile: Interim hail table file name
    :param str truth_htfile: Truth hail table file
    :param int n_bins: Number of bins to create
    :param dict config: Dictionary with threshold values
    :return: Table with bins added
    """
    conf = config["step3"]["create_binned_data_initial"]
    # Count variants for ranking
    # count_expr = {x: hl.agg.filter(hl.is_defined(ht[x]), hl.agg.counter(hl.cond(hl.is_snp(
    #     ht.alleles[0], ht.alleles[1]), 'snv', 'indel'))) for x in ht.row if x.endswith('rank')}
    count_expr = {
        x: hl.agg.filter(
            hl.is_defined(ht[x]), hl.agg.counter(hl.if_else(hl.is_snp(ht.alleles[0], ht.alleles[1]), "snv", "indel"))
        )
        for x in ht.row
        if x.endswith("rank")
    }
    rank_variant_counts = ht.aggregate(hl.Struct(**count_expr))
    print(f"Found the following variant counts:\n {pformat(rank_variant_counts)}")

    ht_truth_data = hl.read_table(path_spark(truth_htfile))
    ht = ht.annotate_globals(rank_variant_counts=fill_empty_rank_values(rank_variant_counts))
    # ht = ht.annotate(
    #     **ht_truth_data[ht.key],
    #     # **fam_ht[ht.key],
    #     # **gnomad_ht[ht.key],
    #     # **denovo_ht[ht.key],
    #     # clinvar=hl.is_defined(clinvar_ht[ht.key]),
    #     indel_length=hl.abs(ht.alleles[0].length()-ht.alleles[1].length()),
    #     rank_bins=hl.array(
    #         [hl.Struct(
    #             rank_id=rank_name,
    #             bin=hl.int(hl.ceil(hl.float(ht[rank_name] + 1) / hl.floor(ht.globals.rank_variant_counts[rank_name][hl.cond(
    #                 hl.is_snp(ht.alleles[0], ht.alleles[1]), 'snv', 'indel')] / n_bins)))
    #         )
    #             for rank_name in rank_variant_counts]
    #     ),
    #     # lcr=hl.is_defined(lcr_intervals[ht.locus])
    # )
    ht = ht.annotate(
        **ht_truth_data[ht.key],
        # **fam_ht[ht.key],
        # **gnomad_ht[ht.key],
        # **denovo_ht[ht.key],
        # clinvar=hl.is_defined(clinvar_ht[ht.key]),
        indel_length=hl.abs(ht.alleles[0].length() - ht.alleles[1].length()),
        rank_bins=hl.array(
            [
                hl.Struct(
                    rank_id=rank_name,
                    bin=hl.int(
                        hl.ceil(
                            hl.float(ht[rank_name] + 1)
                            / hl.floor(
                                ht.globals.rank_variant_counts[rank_name][
                                    hl.if_else(hl.is_snp(ht.alleles[0], ht.alleles[1]), "snv", "indel")
                                ]
                                / n_bins
                            )
                        )
                    ),
                )
                for rank_name in rank_variant_counts
            ]
        ),
        # lcr=hl.is_defined(lcr_intervals[ht.locus])
    )

    ht = ht.explode(ht.rank_bins)
    ht = ht.transmute(rank_id=ht.rank_bins.rank_id, bin=ht.rank_bins.bin)
    ht = ht.filter(hl.is_defined(ht.bin))

    ht = ht.checkpoint(bin_tmp_htfile, overwrite=True)

    # Create binned data
    return (
        ht.group_by(
            rank_id=ht.rank_id,
            contig=ht.locus.contig,
            snv=hl.is_snp(ht.alleles[0], ht.alleles[1]),
            bi_allelic=hl.is_defined(ht.biallelic_rank),
            singleton=ht.transmitted_singleton,
            trans_singletons=hl.is_defined(ht.singleton_rank),
            de_novo_high_quality=ht.de_novo_high_quality_rank,
            de_novo_medium_quality=hl.is_defined(ht.de_novo_medium_quality_rank),
            de_novo_synonymous=hl.is_defined(ht.de_novo_synonymous_rank),
            # release_adj=ht.ac > 0,
            bin=ht.bin,
        )
        ._set_buffer_size(20000)
        .aggregate(
            min_score=hl.agg.min(ht.score),
            max_score=hl.agg.max(ht.score),
            n=hl.agg.count(),
            n_ins=hl.agg.count_where(hl.is_insertion(ht.alleles[0], ht.alleles[1])),
            n_del=hl.agg.count_where(hl.is_deletion(ht.alleles[0], ht.alleles[1])),
            n_ti=hl.agg.count_where(hl.is_transition(ht.alleles[0], ht.alleles[1])),
            n_tv=hl.agg.count_where(hl.is_transversion(ht.alleles[0], ht.alleles[1])),
            n_1bp_indel=hl.agg.count_where(ht.indel_length == 1),
            n_mod3bp_indel=hl.agg.count_where((ht.indel_length % 3) == 0),
            # n_clinvar=hl.agg.count_where(ht.clinvar),
            n_singleton=hl.agg.count_where(ht.transmitted_singleton),
            n_high_quality_de_novos=hl.agg.count_where(ht.de_novo_data.p_de_novo[0] > conf["high_quality_p_de_novo"]),
            # n_validated_DDD_denovos=hl.agg.count_where(
            #    ht.inheritance.contains("De novo")),
            n_medium_quality_de_novos=hl.agg.count_where(
                ht.de_novo_data.p_de_novo[0] > conf["medium_quality_p_de_novo"]
            ),
            n_high_confidence_de_novos=hl.agg.count_where(ht.de_novo_data.confidence[0] == "HIGH"),
            n_de_novo=hl.agg.filter(
                ht.family_stats.unrelated_qc_callstats.AC[0][1] == 0, hl.agg.sum(ht.family_stats.mendel[0].errors)
            ),
            n_high_quality_de_novos_synonymous=hl.agg.count_where(
                (ht.de_novo_data.p_de_novo[0] > conf["high_quality_p_de_novo"])
                & (ht.consequence == "synonymous_variant")
            ),
            n_trans_singletons_synonymous_algorithm=hl.agg.count_where(ht.variant_transmitted_singletons == 1),
            n_untrans_singletons_synonymous_algorithm=hl.agg.count_where(ht.variant_untransmitted_singletons == 1),
            # validated_de_novos=hl.agg.count_where(ht.validated_denovo_inheritance=="De novo constitutive"),
            # n_de_novo_no_lcr=hl.agg.filter(~ht.lcr & (
            #    ht.family_stats.unrelated_qc_callstats.AC[1] == 0), hl.agg.sum(ht.family_stats.mendel.errors)),
            n_de_novo_sites=hl.agg.filter(
                ht.family_stats.unrelated_qc_callstats.AC[0][1] == 0,
                hl.agg.count_where(ht.family_stats.mendel[0].errors > 0),
            ),
            # n_de_novo_sites_no_lcr=hl.agg.filter(~ht.lcr & (
            #    ht.family_stats.unrelated_qc_callstats.AC[1] == 0), hl.agg.count_where(ht.family_stats.mendel.errors > 0)),
            n_trans_singletons=hl.agg.filter(
                (ht.ac_raw < 3) & (ht.family_stats.unrelated_qc_callstats.AC[0][1] == 1),
                hl.agg.sum(ht.family_stats.tdt[0].t),
            ),
            n_trans_singletons_synonymous_hail=hl.agg.filter(
                (ht.ac_raw < 3)
                & (ht.consequence == "synonymous_variant")
                & (ht.family_stats.unrelated_qc_callstats.AC[0][1] == 1),
                hl.agg.sum(ht.family_stats.tdt[0].t),
            ),
            n_untrans_singletons=hl.agg.filter(
                (ht.ac_raw < 3) & (ht.family_stats.unrelated_qc_callstats.AC[0][1] == 1),
                hl.agg.sum(ht.family_stats.tdt[0].u),
            ),
            n_untrans_singletons_synonymous_hail=hl.agg.filter(
                (ht.ac_raw < 3)
                & (ht.consequence == "synonymous_variant")
                & (ht.family_stats.unrelated_qc_callstats.AC[0][1] == 1),
                hl.agg.sum(ht.family_stats.tdt[0].u),
            ),
            n_train_trans_singletons=hl.agg.count_where(
                (ht.family_stats.unrelated_qc_callstats.AC[0][1] == 1) & (ht.family_stats.tdt[0].t == 1)
            ),
            # transmitted and untransmitted common variants
            n_trans_common=hl.agg.filter(
                ht.gnomad_af >= conf["n_trans_common_gnomad_af"], hl.agg.sum(ht.family_stats.tdt[0].t)
            ),
            n_untrans_common=hl.agg.filter(
                ht.gnomad_af >= conf["n_untrans_common_gnomad_af"], hl.agg.sum(ht.family_stats.tdt[0].u)
            ),
            # transmitted/untransmitted with Hail's TDT test
            n_trans_singletons_synonymous_tdt=hl.agg.count_where(
                (ht.consequence == "synonymous_variant")
                & (ht.family_stats.tdt[0].t == 1)
                & (ht.family_stats.tdt[0].u == 0)
            ),
            n_untrans_singletons_synonymous_tdt=hl.agg.count_where(
                (ht.consequence == "synonymous_variant")
                & (ht.family_stats.tdt[0].t == 0)
                & (ht.family_stats.tdt[0].u == 1)
            ),
            # transmitted and untransmitted synonymous variants where unrelated individual allele count <10
            # n_trans_ac_lt_10=hl.agg.filter((ht.family_stats.unrelated_qc_callstats.AC[0][1] <= 10) & (
            #     ht.consequence == "synonymous_variant"), hl.agg.sum(ht.family_stats.tdt[0].t)),
            # n_untrans_ac_lt_10=hl.agg.filter((ht.family_stats.unrelated_qc_callstats.AC[0][1] <= 10) & (
            #     ht.consequence == "synonymous_variant"), hl.agg.sum(ht.family_stats.tdt[0].u)),
            # #transmitted and untransmitted synonymous variants where unrelated individual allele count <7
            # n_trans_ac_lt_7=hl.agg.filter((ht.family_stats.unrelated_qc_callstats.AC[0][1] <= 7) & (
            #     ht.consequence == "synonymous_variant"), hl.agg.sum(ht.family_stats.tdt[0].t)),
            # n_untrans_ac_lt_7=hl.agg.filter((ht.family_stats.unrelated_qc_callstats.AC[0][1] <= 7) & (
            #     ht.consequence == "synonymous_variant"), hl.agg.sum(ht.family_stats.tdt[0].u)),
            # #transmitted and untransmitted synonymous variants where unrelated individual allele count <5
            # n_trans_ac_lt_5=hl.agg.filter((ht.family_stats.unrelated_qc_callstats.AC[0][1] <= 5) & (
            #     ht.consequence == "synonymous_variant"), hl.agg.sum(ht.family_stats.tdt[0].t)),
            # n_untrans_ac_lt_5=hl.agg.filter((ht.family_stats.unrelated_qc_callstats.AC[0][1] <= 5) & (
            #     ht.consequence == "synonymous_variant"), hl.agg.sum(ht.family_stats.tdt[0].u)),
            # #transmitted and untransmitted synonymous variants where unrelated individual allele count <3
            # n_trans_ac_lt_3=hl.agg.filter((ht.family_stats.unrelated_qc_callstats.AC[0][1] <= 3) & (
            #     ht.consequence == "synonymous_variant"), hl.agg.sum(ht.family_stats.tdt[0].t)),
            # n_untrans_ac_lt_3=hl.agg.filter((ht.family_stats.unrelated_qc_callstats.AC[0][1] <= 3) & (
            #     ht.consequence == "synonymous_variant"), hl.agg.sum(ht.family_stats.tdt[0].u)),
            n_omni=hl.agg.count_where(ht.omni),
            n_mills=hl.agg.count_where(ht.mills),
            n_hapmap=hl.agg.count_where(ht.hapmap),
            n_kgp_high_conf_snvs=hl.agg.count_where(ht.kgp_phase1_hc),
            fail_hard_filters=hl.agg.count_where(ht.fail_hard_filters),
            fail_hard_filters_snvs=hl.agg.count_where(
                (ht.fail_hard_filters) & (hl.is_snp(ht.alleles[0], ht.alleles[1]))
            ),
            fail_hard_filters_indels=hl.agg.count_where(
                (ht.fail_hard_filters) & (hl.is_indel(ht.alleles[0], ht.alleles[1]))
            ),
            # n_vqsr_pos_train=hl.agg.count_where(ht.vqsr_positive_train_site),
            # n_vqsr_neg_train=hl.agg.count_where(ht.vqsr_negative_train_site)
        )
    )


def fill_empty_rank_values(rank_counts: hl.Struct) -> hl.Struct:
    filled_rank_counts = {}
    for field in rank_counts:
        if rank_counts[field] == {}:
            filled_rank_counts[field] = {"indel": 0, "snv": 0}
        else:
            filled_rank_counts[field] = rank_counts[field]
    return hl.struct(**filled_rank_counts)


def make_subrank_expression(ht: hl.Table, add_trios_dependent_subranks: bool) -> Dict[str, hl.expr.BooleanExpression]:
    # TODO: Do we need to move hardcoded values to config?

    subrank_expr = {
        "singleton_rank": ht.transmitted_singleton,
        "biallelic_rank": ~ht.was_split,
        "biallelic_singleton_rank": ~ht.was_split & ht.transmitted_singleton,
        "de_novo_high_quality_rank": ht.de_novo_data.p_de_novo[0] > 0.9,
        "de_novo_medium_quality_rank": ht.de_novo_data.p_de_novo[0] > 0.5,
        "de_novo_synonymous_rank": ht.consequence == "synonymous_variant",
    }

    if not add_trios_dependent_subranks:
        # We know that there will be no de-novo and singleton annotations: change all fields to False
        subrank_expr["singleton_rank"] = hl.literal(False)
        subrank_expr["biallelic_singleton_rank"] = hl.literal(False)
        subrank_expr["de_novo_high_quality_rank"] = hl.literal(False)
        subrank_expr["de_novo_medium_quality_rank"] = hl.literal(False)
    return subrank_expr


def main():
    # = STEP SETUP =
    config = parse_config()
    tmp_dir = config["general"]["tmp_dir"]

    # = STEP PARAMETERS = #
    model_id = config["general"]["rf_model_id"]
    pedfile: str = config["step3"]["pedfile"]

    # = STEP DEPENDENCIES = #
    rf_dir: str = path_spark(config["general"]["var_qc_rf_dir"])
    # TODO: should be in the config
    htfile: str = os.path.join(rf_dir, model_id, "rf_result_final_for_ranking.ht")

    # = STEP OUTPUTS = #
    htrankedfile: str = os.path.join(rf_dir, model_id, "rf_result_ranked.ht")
    truth_htfile: str = config["step3"]["create_binned_data_initial"]["truth_htfile"]
    # TODO: should be in the config
    bin_tmp_htfile: str = os.path.join(rf_dir, model_id, "_gnomad_score_binning_tmp.ht")
    bin_htfile: str = os.path.join(rf_dir, model_id, "_rf_result_ranked_BINS.ht")

    # = STEP LOGIC = #
    _ = hail_utils.init_hl(tmp_dir)
    ht = hl.read_table(path_spark(htfile))
    print("=== Assigning ranks ===")

    add_trios_dependent_subranks = pedfile is not None
    subrank_expr = make_subrank_expression(ht, add_trios_dependent_subranks)

    ht_ranked = add_rank(
        ht,
        score_expr=(1 - ht.rf_probability["TP"]),
        subrank_expr=subrank_expr,
    )
    ht_ranked = ht_ranked.annotate(score=(1 - ht_ranked.rf_probability["TP"]))
    ht_ranked.write(path_spark(htrankedfile), overwrite=True)

    print("=== Ranking results and adding bins ===")
    ht_bins = create_binned_data_initial(ht_ranked, bin_tmp_htfile, truth_htfile, n_bins=100, config=config)
    ht_bins.write(path_spark(bin_htfile), overwrite=True)


if __name__ == "__main__":
    main()
