from typing import Tuple, Optional

import hail as hl
import os.path

from wes_qc import hail_utils
from utils.utils import parse_config, path_spark


def add_cq_annotation(ht: hl.Table, synonymous_file: str) -> hl.Table:
    """
    Add annotation for synonymous variants
    ;param str htfile: Random forest hail table file
    ;param str synonymous_file: text file containing synonymous variants
    ;param str ht_cq_fle: Output hail table file annotated with synonymous variants
    """
    # DEBUG: test if file:// prefix is needed for this function
    synonymous_ht = hl.import_table(
        path_spark(synonymous_file),
        types={"f0": "str", "f1": "int32", "f2": "str", "f3": "str", "f4": "str", "f5": "str"},
        no_header=True,
    )
    synonymous_ht = synonymous_ht.annotate(chr=synonymous_ht.f0)
    synonymous_ht = synonymous_ht.annotate(pos=synonymous_ht.f1)
    synonymous_ht = synonymous_ht.annotate(rs=synonymous_ht.f2)
    synonymous_ht = synonymous_ht.annotate(ref=synonymous_ht.f3)
    synonymous_ht = synonymous_ht.annotate(alt=synonymous_ht.f4)
    synonymous_ht = synonymous_ht.annotate(consequence=synonymous_ht.f5)
    synonymous_ht = synonymous_ht.key_by(
        locus=hl.locus(synonymous_ht.chr, synonymous_ht.pos), alleles=[synonymous_ht.ref, synonymous_ht.alt]
    )
    synonymous_ht = synonymous_ht.drop(
        synonymous_ht.f0,
        synonymous_ht.f1,
        synonymous_ht.f2,
        synonymous_ht.f3,
        synonymous_ht.f4,
        synonymous_ht.chr,
        synonymous_ht.pos,
        synonymous_ht.ref,
        synonymous_ht.alt,
    )
    synonymous_ht = synonymous_ht.key_by(synonymous_ht.locus, synonymous_ht.alleles)

    ht_cq = ht.annotate(consequence=synonymous_ht[ht.key].consequence)
    return ht_cq


def dnm_and_family_annotation_with_missing(ht_cq: hl.Table) -> hl.Table:
    """
    Annotates the input Hail Table with missing family annotations and returns the updated table.
    """
    # Hail types for missing family annotations
    # You can obtain it from any real family stats table using the .dtype property of the column
    de_novo_data_dtype = "array<struct{id: str, prior: float64, proband: struct{s: str, batch: str, c_a_artefact_suspected: bool, assigned_pop: str}, father: struct{s: str, batch: str, c_a_artefact_suspected: bool, assigned_pop: str}, mother: struct{s: str, batch: str, c_a_artefact_suspected: bool, assigned_pop: str}, proband_entry: struct{AD: array<int32>, DP: int32, GQ: int32, GT: call, MIN_DP: int32, PGT: call, PID: str, PL: array<int32>, PS: int32, RGQ: int32, SB: array<int32>, sum_AD: int32, adj: bool, HetAB: float64}, father_entry: struct{AD: array<int32>, DP: int32, GQ: int32, GT: call, MIN_DP: int32, PGT: call, PID: str, PL: array<int32>, PS: int32, RGQ: int32, SB: array<int32>, sum_AD: int32, adj: bool, HetAB: float64}, mother_entry: struct{AD: array<int32>, DP: int32, GQ: int32, GT: call, MIN_DP: int32, PGT: call, PID: str, PL: array<int32>, PS: int32, RGQ: int32, SB: array<int32>, sum_AD: int32, adj: bool, HetAB: float64}, is_female: bool, p_de_novo: float64, confidence: str}>"
    family_stats_dtype = "array<struct{mendel: struct{errors: int64}, tdt: struct{t: int64, u: int64, chi_sq: float64, p_value: float64}, unrelated_qc_callstats: struct{AC: array<int32>, AF: array<float64>, AN: int32, homozygote_count: array<int32>}, meta: dict<str, str>}>"
    trio_stat_dtype = "struct{n_transmitted_raw: int64, n_untransmitted_raw: int64, n_transmitted_adj: int64, n_untransmitted_adj: int64, n_de_novos_raw: int64, n_de_novos_adj: int64, ac_parents_raw: int64, an_parents_raw: int64, ac_children_raw: int64, an_children_raw: int64, ac_parents_adj: int64, an_parents_adj: int64, ac_children_adj: int64, an_children_adj: int64}"

    ht_cq = ht_cq.annotate(de_novo_data=hl.missing(de_novo_data_dtype))
    ht_cq = ht_cq.annotate(family_stats=hl.missing(family_stats_dtype))
    ht_cq = ht_cq.annotate(fam=hl.missing(trio_stat_dtype))

    return ht_cq


def dnm_and_family_annotation_with_tables(
    ht_cq: hl.Table, dnm_htfile: str, fam_stats_htfile: str, trio_stats_htfile: str
):
    """
    Annotate RF result hail table with DNMs and family stats
    """

    # Reading annotations or mocking it
    de_novo_annot = hl.read_table(path_spark(dnm_htfile))[ht_cq.key].de_novo_data
    ht_cq = ht_cq.annotate(de_novo_data=de_novo_annot)

    family_stat_annot = hl.read_table(path_spark(fam_stats_htfile))[ht_cq.key].family_stats
    ht_cq = ht_cq.annotate(family_stats=family_stat_annot)

    trio_stat_annot = hl.read_table(path_spark(trio_stats_htfile))[ht_cq.key]
    ht_cq = ht_cq.annotate(fam=trio_stat_annot)

    return ht_cq


def dnm_and_family_annotation(
    ht_cq: hl.Table, dnm_htfile: str, fam_stats_htfile: str, trio_stats_htfile: str, pedfile: str
) -> hl.Table:
    if pedfile is not None:
        # annotate with family stats and DNMs
        ht_cq = dnm_and_family_annotation_with_tables(ht_cq, dnm_htfile, fam_stats_htfile, trio_stats_htfile)
    else:
        ht_cq = dnm_and_family_annotation_with_missing(ht_cq)
    return ht_cq


def count_trans_untransmitted_singletons(mt_trios: Optional[hl.MatrixTable], ht: hl.Table) -> Tuple[int, int]:
    """
    Count untransmitted and transmitted singletons
    """
    if mt_trios is not None:
        # Filtering trios matrixtable
        mt_trios = mt_trios.annotate_rows(consequence=ht[mt_trios.row_key].consequence)

        # TODO: there is a save step here in Pavlos file, is it needed?
        # mt_filtered = mt_trios.filter_rows((mt_trios.variant_qc.AC[1] <= 2) & (mt_trios.consequence == "synonymous_variant"))
        # mt_filtered = mt_trios.filter_entries((mt_trios.variant_qc.AC[1] <= 2) & (mt_trios.consequence == "synonymous_variant"))
        # TODO: unused step - was not commented out - probably incorrect
        # mt_filtered = mt_trios.filter_rows(
        #     (mt_trios.varqc_trios.AC[1] <= 2) & (mt_trios.consequence == "synonymous_variant")
        # )

        mt_filtered = mt_trios.filter_entries(
            (mt_trios.varqc_trios.AC[1] <= 2) & (mt_trios.consequence == "synonymous_variant")
        )

        mt_trans = mt_filtered.filter_entries(mt_filtered.varqc_trios.AC[1] == 2)
        mt_untrans = mt_filtered.filter_entries(mt_filtered.varqc_trios.AC[1] == 1)

        mt_trans_count = mt_trans.group_cols_by(mt_trans.id).aggregate(
            transmitted_singletons_count=hl.agg.count_where(
                # (mt_trans.info.AC[0] == 2) &
                (mt_trans.proband_entry.GT.is_non_ref())
                & ((mt_trans.father_entry.GT.is_non_ref()) | (mt_trans.mother_entry.GT.is_non_ref()))
            )
        )

        total_transmitted_singletons = mt_trans_count.aggregate_entries(
            hl.agg.count_where(mt_trans_count.transmitted_singletons_count > 0)
        )

        mt_untrans_count = mt_untrans.group_cols_by(mt_untrans.id).aggregate(
            untransmitted_singletons_count=hl.agg.count_where(
                # (mt_untrans.info.AC[0] == 1) &
                (mt_untrans.proband_entry.GT.is_hom_ref())
                & ((mt_untrans.father_entry.GT.is_non_ref()) | (mt_untrans.mother_entry.GT.is_non_ref()))
            )
        )
        total_untransmitted_singletons = mt_untrans_count.aggregate_entries(
            hl.agg.count_where(mt_untrans_count.untransmitted_singletons_count > 0)
        )

        print(f"\nTransmitted singletons:{total_transmitted_singletons}\n")
        print(f"\nUntransmitted singletons:{total_untransmitted_singletons}\n")

        if total_untransmitted_singletons > 0:
            ratio_transmitted_untransmitted = total_transmitted_singletons / total_untransmitted_singletons
            print(f"\nTrans/Untrans ratio:{ratio_transmitted_untransmitted}\n")
        mt2 = mt_trans_count.annotate_rows(
            variant_transmitted_singletons=hl.agg.count_where(mt_trans_count.transmitted_singletons_count == 1)
        )
        mt2.variant_transmitted_singletons.summarize()

        mt3 = mt_untrans_count.annotate_rows(
            variant_untransmitted_singletons=hl.agg.count_where(mt_untrans_count.untransmitted_singletons_count == 1)
        )
        mt3.variant_untransmitted_singletons.summarize()

        variant_transmitted_annot = mt2.rows()[ht.key].variant_transmitted_singletons
        ht = ht.annotate(variant_transmitted_singletons=variant_transmitted_annot)

        variant_untransmitted_annot = mt3.rows()[ht.key].variant_untransmitted_singletons
        ht = ht.annotate(variant_untransmitted_singletons=variant_untransmitted_annot)
    else:
        print("=== No pedifree data found: skipping transmitted/untransmitted calculations")
        ht = ht.annotate(variant_transmitted_singletons=0, variant_untransmitted_singletons=0)
    return ht


def transmitted_singleton_annotation(
    ht: hl.Table,
    trio_mtfile: str,
    pedfile: str,
) -> hl.Table:
    """
    Annotate MT with transmited singletons and create final variant QC HT for ranking
    """
    print("Annotating with transmitted singletons")
    mt_trios = hl.read_matrix_table(path_spark(trio_mtfile)) if pedfile is not None else None
    ht = count_trans_untransmitted_singletons(mt_trios, ht)
    return ht


def annotate_gnomad_af(ht: hl.Table, gnomad_ht: hl.Table) -> hl.Table:
    """
    Annotate with gnomad allele frequencies
    """
    print("=== Annotating with gnomad AF ===")
    ht = ht.annotate(gnomad_af=gnomad_ht[ht.key].freq[0].AF)
    return ht


def main():
    # = STEP SETUP = #
    config = parse_config()
    tmp_dir = config["general"]["tmp_dir"]

    # = STEP PARAMETERS = #
    pedfile: str = config["step3"]["pedfile"]
    model_id: str = config["general"]["rf_model_id"]
    synonymous_file: str = config["step3"]["add_cq_annotation"]["synonymous_file"]
    gnomad_htfile: str = config["step3"]["annotate_gnomad"]["gnomad_htfile"]

    # = STEP DEPENDENCIES = #
    rf_dir = path_spark(config["general"]["var_qc_rf_dir"])
    htfile: str = os.path.join(rf_dir, model_id, "rf_result.ht")

    # = STEP OUTPUTS = #
    dnm_htfile: str = config["step3"]["dnm_and_family_annotation"]["dnm_htfile"]
    fam_stats_htfile: str = config["step3"]["dnm_and_family_annotation"]["fam_stats_htfile"]
    trio_stats_htfile: str = config["step3"]["dnm_and_family_annotation"]["trio_stats_htfile"]
    family_annot_htfile: str = os.path.join(rf_dir, model_id, "rf_result_denovo_family_stats.ht")
    trio_mtfile: str = config["step3"]["transmitted_singleton_annotation"]["trio_mtfile"]
    # trio_filtered_mtfile:str = config["step3"]["transmitted_singleton_annotation"]["trio_filtered_mtfile"] # TODO: remove form config if not needed
    trans_sing_htfile: str = os.path.join(rf_dir, model_id, "rf_result_trans_sing.ht")
    ht_cq_file: str = os.path.join(rf_dir, model_id, "rf_result_with_synonymous.ht")
    final_htfile: str = os.path.join(rf_dir, model_id, "rf_result_final_for_ranking.ht")

    # = STEP LOGIC = #
    _ = hail_utils.init_hl(tmp_dir)

    # annotate with synonymous CQs
    ht = hl.read_table(path_spark(htfile))
    ht_cq = add_cq_annotation(ht, synonymous_file)
    ht_cq.write(path_spark(ht_cq_file), overwrite=True)

    # ht_cq = hl.read_table(path_spark(ht_cq_file))
    ht_cq = dnm_and_family_annotation(ht_cq, dnm_htfile, fam_stats_htfile, trio_stats_htfile, pedfile)
    ht_cq.write(path_spark(family_annot_htfile), overwrite=True)

    # annotate with transmitted singletons
    # ht = hl.read_table(path_spark(family_annot_htfile))
    ht_cq = transmitted_singleton_annotation(ht_cq, trio_mtfile, pedfile)
    ht_cq.write(path_spark(trans_sing_htfile), overwrite=True)

    # annotate with gnomad AF
    gnomad_ht = hl.read_table(path_spark(gnomad_htfile))
    ht_cq = annotate_gnomad_af(ht_cq, gnomad_ht)
    ht_cq.write(final_htfile, overwrite=True)


if __name__ == "__main__":
    main()
