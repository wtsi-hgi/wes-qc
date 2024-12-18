import hail as hl
import os.path
from wes_qc import hail_utils
from utils.utils import parse_config, rm_mt, path_spark


def add_cq_annotation(htfile: str, synonymous_file: str, ht_cq_file: str):
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

    ht = hl.read_table(path_spark(htfile))
    ht = ht.annotate(consequence=synonymous_ht[ht.key].consequence)
    ht.write(path_spark(ht_cq_file), overwrite=True)


def dnm_and_family_annotation(
    ht_cq_fle: str, dnm_htfile: str, fam_stats_htfile: str, trio_stats_htfile: str, family_annot_htfile: str
):
    """
    Annotate RF result hail table with DNMs and family stats
    ;param str ht_cq_fle: Input hail table file annotated with synonymous variants
    :param str dnm_htfile: De novo hail table file
    :param str fam_stats_htfile: Family stats hail table file
    :param str trio_stats_htfile: Trio stats hail table file
    :param str family_annot_htfile: Output hail table annotated with DNMs and family stats file
    """
    ht = hl.read_table(path_spark(ht_cq_fle))
    dnm_ht = hl.read_table(path_spark(dnm_htfile))
    fam_stats_ht = hl.read_table(path_spark(fam_stats_htfile))
    trio_stats_ht = hl.read_table(path_spark(trio_stats_htfile))
    ht = ht.annotate(de_novo_data=dnm_ht[ht.key].de_novo_data)
    ht = ht.annotate(family_stats=fam_stats_ht[ht.key].family_stats)
    ht = ht.annotate(fam=trio_stats_ht[ht.key])
    ht.write(path_spark(family_annot_htfile), overwrite=True)


def count_trans_untransmitted_singletons(mt_filtered: hl.MatrixTable, ht: hl.Table) -> hl.Table:
    """
    Count untransmitted and transmitted singletons
    :param hl.MatrixTable mt_filtered: Filtered trio matrixtable
    :param hl.Table ht: Output hail table
    """

    # mt_trans = mt_filtered.filter_entries(mt_filtered.variant_qc.AC[1] == 2)
    # mt_untrans = mt_filtered.filter_entries(mt_filtered.variant_qc.AC[1] == 1)
    mt_trans = mt_filtered.filter_entries(mt_filtered.varqc_trios.AC[1] == 2)
    mt_untrans = mt_filtered.filter_entries(mt_filtered.varqc_trios.AC[1] == 1)

    mt_trans_count = mt_trans.group_cols_by(mt_trans.id).aggregate(
        transmitted_singletons_count=hl.agg.count_where(
            # (mt_trans.info.AC[0] == 2) &
            (mt_trans.proband_entry.GT.is_non_ref())
            & ((mt_trans.father_entry.GT.is_non_ref()) | (mt_trans.mother_entry.GT.is_non_ref()))
        )
    )

    Total_transmitted_singletons = mt_trans_count.aggregate_entries(
        hl.agg.count_where(mt_trans_count.transmitted_singletons_count > 0)
    )
    print(Total_transmitted_singletons)

    mt_untrans_count = mt_untrans.group_cols_by(mt_untrans.id).aggregate(
        untransmitted_singletons_count=hl.agg.count_where(
            # (mt_untrans.info.AC[0] == 1) &
            (mt_untrans.proband_entry.GT.is_hom_ref())
            & ((mt_untrans.father_entry.GT.is_non_ref()) | (mt_untrans.mother_entry.GT.is_non_ref()))
        )
    )
    Total_untransmitted_singletons = mt_untrans_count.aggregate_entries(
        hl.agg.count_where(mt_untrans_count.untransmitted_singletons_count > 0)
    )
    print(Total_untransmitted_singletons)

    print(f"\nTransmitted singletons:{Total_transmitted_singletons}\n")
    print(f"\nUntransmitted singletons:{Total_untransmitted_singletons}")

    if Total_untransmitted_singletons > 0:
        Ratio_transmitted_untransmitted = Total_transmitted_singletons / Total_untransmitted_singletons
        print(f"\nRatio:{Ratio_transmitted_untransmitted}\n")
    mt2 = mt_trans_count.annotate_rows(
        variant_transmitted_singletons=hl.agg.count_where(mt_trans_count.transmitted_singletons_count == 1)
    )
    mt2.variant_transmitted_singletons.summarize()

    mt3 = mt_untrans_count.annotate_rows(
        variant_untransmitted_singletons=hl.agg.count_where(mt_untrans_count.untransmitted_singletons_count == 1)
    )
    mt3.variant_untransmitted_singletons.summarize()

    ht = ht.annotate(variant_transmitted_singletons=mt2.rows()[ht.key].variant_transmitted_singletons)
    ht = ht.annotate(variant_untransmitted_singletons=mt3.rows()[ht.key].variant_untransmitted_singletons)
    return ht


def transmitted_singleton_annotation(
    family_annot_htfile: str, trio_mtfile: str, trio_filtered_mtfile: str, trans_sing_htfile: str
):
    """
    Annotate MT with transmited singletons and create final variant QC HT for ranking
    :param str family_annot_htfile: Family annotation hail table file
    :param str trio_mtfile: Trio annotation hail matrixtable file
    :param str trio_filtered_mtfile: Trio filtered hail matrixtable file
    :param str trans_sing_htfile: Variant QC hail table file with transmitted singleton annotation
    """
    print("Annotating with transmitted singletons")
    ht = hl.read_table(path_spark(family_annot_htfile))
    mt_trios = hl.read_matrix_table(path_spark(trio_mtfile))
    mt_trios = mt_trios.annotate_rows(consequence=ht[mt_trios.row_key].consequence)
    # there is a save step here in Pavlos file, is it needed?
    # mt_filtered = mt_trios.filter_rows((mt_trios.variant_qc.AC[1] <= 2) & (mt_trios.consequence == "synonymous_variant"))
    # mt_filtered = mt_trios.filter_entries((mt_trios.variant_qc.AC[1] <= 2) & (mt_trios.consequence == "synonymous_variant"))
    mt_filtered = mt_trios.filter_rows(
        (mt_trios.varqc_trios.AC[1] <= 2) & (mt_trios.consequence == "synonymous_variant")
    )
    mt_filtered = mt_trios.filter_entries(
        (mt_trios.varqc_trios.AC[1] <= 2) & (mt_trios.consequence == "synonymous_variant")
    )
    mt_filtered = mt_filtered.checkpoint(path_spark(trio_filtered_mtfile), overwrite=True)

    ht = count_trans_untransmitted_singletons(mt_filtered, ht)
    ht.write(path_spark(trans_sing_htfile), overwrite=True)
    rm_mt(trio_filtered_mtfile)


def annotate_gnomad(tdt_htfile: str, gnomad_htfile: str, final_htfile: str):
    """
    Annotate with gnomad allele frequencies
    :param str tdt_htfile: Htfile with transmitted/untransmitted counts
    :param str gnomad_htfile: Gnomad annotation hail table file
    :param str final_htfile: Final RF htfile for ranking and binning
    """
    print("Annotating with gnomad AF")
    ht = hl.read_table(path_spark(tdt_htfile))
    gnomad_ht = hl.read_table(path_spark(gnomad_htfile))
    ht = ht.annotate(gnomad_af=gnomad_ht[ht.key].freq[0].AF)
    ht.write(final_htfile, overwrite=True)


def main():
    # = STEP SETUP = #
    config = parse_config()
    tmp_dir = config["general"]["tmp_dir"]

    # = STEP PARAMETERS = #
    model_id = config["general"]["rf_model_id"]
    synonymous_file = config["step3"]["add_cq_annotation"]["synonymous_file"]
    gnomad_htfile = config["step3"]["annotate_gnomad"]["gnomad_htfile"]

    # = STEP DEPENDENCIES = #
    rf_dir = path_spark(config["general"]["var_qc_rf_dir"])
    # TODO: should be in the config
    htfile = os.path.join(rf_dir, model_id, "rf_result.ht")

    # = STEP OUTPUTS = #
    dnm_htfile = config["step3"]["dnm_and_family_annotation"]["dnm_htfile"]
    fam_stats_htfile = config["step3"]["dnm_and_family_annotation"]["fam_stats_htfile"]
    trio_stats_htfile = config["step3"]["dnm_and_family_annotation"]["trio_stats_htfile"]
    family_annot_htfile = os.path.join(rf_dir, model_id, "rf_result_denovo_family_stats.ht")
    trio_mtfile = config["step3"]["transmitted_singleton_annotation"]["trio_mtfile"]
    trio_filtered_mtfile = config["step3"]["transmitted_singleton_annotation"]["trio_filtered_mtfile"]
    trans_sing_htfile = os.path.join(rf_dir, model_id, "rf_result_trans_sing.ht")
    # TODO: should be in the config
    ht_cq_file = os.path.join(rf_dir, model_id, "rf_result_with_synonymous.ht")
    # TODO: should be in the config
    final_htfile = os.path.join(rf_dir, model_id, "rf_result_final_for_ranking.ht")

    # = STEP LOGIC = #
    _ = hail_utils.init_hl(tmp_dir)

    # annotate with synonymous CQs
    add_cq_annotation(htfile, synonymous_file, ht_cq_file)

    # annotate with family stats and DNMs
    dnm_and_family_annotation(ht_cq_file, dnm_htfile, fam_stats_htfile, trio_stats_htfile, family_annot_htfile)

    # annotate with transmitted singletons
    transmitted_singleton_annotation(family_annot_htfile, trio_mtfile, trio_filtered_mtfile, trans_sing_htfile)

    # annotate with gnomad AF
    annotate_gnomad(trans_sing_htfile, gnomad_htfile, final_htfile)


if __name__ == "__main__":
    main()
