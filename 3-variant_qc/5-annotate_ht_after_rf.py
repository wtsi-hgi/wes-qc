from typing import Optional

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
    # Define Hail type objects for missing family annotations
    de_novo_data_type = hl.tarray(
        hl.tstruct(
            id=hl.tstr,
            prior=hl.tfloat64,
            proband=hl.tstruct(s=hl.tstr, batch=hl.tstr, c_a_artefact_suspected=hl.tbool, assigned_pop=hl.tstr),
            father=hl.tstruct(s=hl.tstr, batch=hl.tstr, c_a_artefact_suspected=hl.tbool, assigned_pop=hl.tstr),
            mother=hl.tstruct(s=hl.tstr, batch=hl.tstr, c_a_artefact_suspected=hl.tbool, assigned_pop=hl.tstr),
            proband_entry=hl.tstruct(
                AD=hl.tarray(hl.tint32),
                DP=hl.tint32,
                GQ=hl.tint32,
                GT=hl.tcall,
                MIN_DP=hl.tint32,
                PGT=hl.tcall,
                PID=hl.tstr,
                PL=hl.tarray(hl.tint32),
                PS=hl.tint32,
                RGQ=hl.tint32,
                SB=hl.tarray(hl.tint32),
                sum_AD=hl.tint32,
                adj=hl.tbool,
                HetAB=hl.tfloat64,
            ),
            father_entry=hl.tstruct(
                AD=hl.tarray(hl.tint32),
                DP=hl.tint32,
                GQ=hl.tint32,
                GT=hl.tcall,
                MIN_DP=hl.tint32,
                PGT=hl.tcall,
                PID=hl.tstr,
                PL=hl.tarray(hl.tint32),
                PS=hl.tint32,
                RGQ=hl.tint32,
                SB=hl.tarray(hl.tint32),
                sum_AD=hl.tint32,
                adj=hl.tbool,
                HetAB=hl.tfloat64,
            ),
            mother_entry=hl.tstruct(
                AD=hl.tarray(hl.tint32),
                DP=hl.tint32,
                GQ=hl.tint32,
                GT=hl.tcall,
                MIN_DP=hl.tint32,
                PGT=hl.tcall,
                PID=hl.tstr,
                PL=hl.tarray(hl.tint32),
                PS=hl.tint32,
                RGQ=hl.tint32,
                SB=hl.tarray(hl.tint32),
                sum_AD=hl.tint32,
                adj=hl.tbool,
                HetAB=hl.tfloat64,
            ),
            is_female=hl.tbool,
            p_de_novo=hl.tfloat64,
            confidence=hl.tstr,
        )
    )

    family_stats_type = hl.tarray(
        hl.tstruct(
            mendel=hl.tstruct(errors=hl.tint64),
            tdt=hl.tstruct(t=hl.tint64, u=hl.tint64, chi_sq=hl.tfloat64, p_value=hl.tfloat64),
            unrelated_qc_callstats=hl.tstruct(
                AC=hl.tarray(hl.tint32), AF=hl.tarray(hl.tfloat64), AN=hl.tint32, homozygote_count=hl.tarray(hl.tint32)
            ),
            meta=hl.tdict(hl.tstr, hl.tstr),
        )
    )

    trio_stat_type = hl.tstruct(
        n_transmitted_raw=hl.tint64,
        n_untransmitted_raw=hl.tint64,
        n_transmitted_adj=hl.tint64,
        n_untransmitted_adj=hl.tint64,
        n_de_novos_raw=hl.tint64,
        n_de_novos_adj=hl.tint64,
        ac_parents_raw=hl.tint64,
        an_parents_raw=hl.tint64,
        ac_children_raw=hl.tint64,
        an_children_raw=hl.tint64,
        ac_parents_adj=hl.tint64,
        an_parents_adj=hl.tint64,
        ac_children_adj=hl.tint64,
        an_children_adj=hl.tint64,
    )

    ht_cq = ht_cq.annotate(de_novo_data=hl.missing(de_novo_data_type))
    ht_cq = ht_cq.annotate(family_stats=hl.missing(family_stats_type))
    ht_cq = ht_cq.annotate(fam=hl.missing(trio_stat_type))

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
    ht: hl.Table, dnm_htfile: str, fam_stats_htfile: str, trio_stats_htfile: str, pedfile: str
) -> hl.Table:
    if pedfile is not None:
        # annotate with family stats and DNMs
        print("=== Annotating with family stats and DNMs ===")
        ht = dnm_and_family_annotation_with_tables(ht, dnm_htfile, fam_stats_htfile, trio_stats_htfile)
    else:
        print("=== No pedifree data found: annotating with missing family stats and DNMs ===")
        ht = dnm_and_family_annotation_with_missing(ht)
    return ht


# TODO: To messy - a good candidate for refactoring
def count_trans_untransmitted_singletons(mt_trios: Optional[hl.MatrixTable], ht: hl.Table) -> hl.Table:
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
    ht = add_cq_annotation(ht, synonymous_file)
    ht.write(path_spark(ht_cq_file), overwrite=True)

    # ht = hl.read_table(path_spark(ht_cq_file))
    ht = dnm_and_family_annotation(ht, dnm_htfile, fam_stats_htfile, trio_stats_htfile, pedfile)
    ht.write(path_spark(family_annot_htfile), overwrite=True)

    # annotate with transmitted singletons
    # ht = hl.read_table(path_spark(family_annot_htfile))
    ht = transmitted_singleton_annotation(ht, trio_mtfile, pedfile)
    ht.write(path_spark(trans_sing_htfile), overwrite=True)

    # annotate with gnomad AF
    gnomad_ht = hl.read_table(path_spark(gnomad_htfile))
    ht = annotate_gnomad_af(ht, gnomad_ht)
    ht.write(final_htfile, overwrite=True)


if __name__ == "__main__":
    main()
