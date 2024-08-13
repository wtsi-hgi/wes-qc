#
import os
from typing import Any, Optional

import hail as hl
import pyspark
import argparse
from utils.utils import parse_config
from wes_qc import hail_utils


def add_cq_annotation(htfile: str, synonymous_file: str, ht_cq_file: str) -> None:
    """
    Add annotation for synonymous variants
    ;param str htfile: Random forest hail table file
    ;param str synonymous_file: text file containing synonymous variants
    ;param str ht_cq_fle: Output hail table file annotated with synonymous variants
    """
    synonymous_ht = hl.import_table(
        synonymous_file,
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

    ht = hl.read_table(htfile)
    ht = ht.annotate(consequence=synonymous_ht[ht.key].consequence)
    ht.write(ht_cq_file, overwrite=True)


def dnm_and_family_annotation(
    ht_cq_fle: str,
    dnm_htfile: str,
    fam_stats_htfile: str,
    trio_stats_htfile: str,
    family_annot_htfile: str,
    mock: bool = False,
) -> None:
    """
    Annotate RF result hail table with DNMs and family stats
    ;param str ht_cq_fle: Input hail table file annotated with synonymous variants
    :param str dnm_htfile: De novo hail table file
    :param str fam_stats_htfile: Family stats hail table file
    :param str trio_stats_htfile: Trio stats hail table file
    :param str family_annot_htfile: Output hail table annotated with DNMs and family stats file
    """

    # The types for missing family annotations
    # You can obtain it from any real family stats table using the .dtype property of the column
    de_novo_data_dtype = "array<struct{id: str, prior: float64, proband: struct{s: str, batch: str, c_a_artefact_suspected: bool, assigned_pop: str}, father: struct{s: str, batch: str, c_a_artefact_suspected: bool, assigned_pop: str}, mother: struct{s: str, batch: str, c_a_artefact_suspected: bool, assigned_pop: str}, proband_entry: struct{AD: array<int32>, DP: int32, GQ: int32, GT: call, MIN_DP: int32, PGT: call, PID: str, PL: array<int32>, PS: int32, RGQ: int32, SB: array<int32>, sum_AD: int32, adj: bool, HetAB: float64}, father_entry: struct{AD: array<int32>, DP: int32, GQ: int32, GT: call, MIN_DP: int32, PGT: call, PID: str, PL: array<int32>, PS: int32, RGQ: int32, SB: array<int32>, sum_AD: int32, adj: bool, HetAB: float64}, mother_entry: struct{AD: array<int32>, DP: int32, GQ: int32, GT: call, MIN_DP: int32, PGT: call, PID: str, PL: array<int32>, PS: int32, RGQ: int32, SB: array<int32>, sum_AD: int32, adj: bool, HetAB: float64}, is_female: bool, p_de_novo: float64, confidence: str}>"
    family_stats_dtype = "array<struct{mendel: struct{errors: int64}, tdt: struct{t: int64, u: int64, chi_sq: float64, p_value: float64}, unrelated_qc_callstats: struct{AC: array<int32>, AF: array<float64>, AN: int32, homozygote_count: array<int32>}, meta: dict<str, str>}>"
    trio_stat_dtype = "struct{n_transmitted_raw: int64, n_untransmitted_raw: int64, n_transmitted_adj: int64, n_untransmitted_adj: int64, n_de_novos_raw: int64, n_de_novos_adj: int64, ac_parents_raw: int64, an_parents_raw: int64, ac_children_raw: int64, an_children_raw: int64, ac_parents_adj: int64, an_parents_adj: int64, ac_children_adj: int64, an_children_adj: int64}"

    ht = hl.read_table(ht_cq_fle)

    de_novo_annot = hl.read_table(dnm_htfile)[ht.key].de_novo_data if not mock else hl.missing(de_novo_data_dtype)
    ht = ht.annotate(de_novo_data=de_novo_annot)

    family_stat_annot = (
        hl.read_table(fam_stats_htfile)[ht.key].family_stats if not mock else hl.missing(family_stats_dtype)
    )
    ht = ht.annotate(family_stats=family_stat_annot)

    trio_stat_annot = hl.read_table(trio_stats_htfile)[ht.key] if not mock else hl.missing(trio_stat_dtype)
    ht = ht.annotate(fam=trio_stat_annot)

    ht.write(family_annot_htfile, overwrite=True)


def count_trans_untransmitted_singletons(mt_filtered: Optional[hl.MatrixTable], ht: hl.Table) -> hl.Table:
    """
    Count untransmitted and transmitted singletons
    :param hl.MatrixTable mt_filtered: Filtered trio matrixtable
    :param hl.Table ht: Output hail table
    """
    if mt_filtered is not None:
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
        print(f"\nUntransmitted singletons:{total_untransmitted_singletons}")

        Ratio_transmitted_untransmitted = total_transmitted_singletons / total_untransmitted_singletons
        print(Ratio_transmitted_untransmitted)
        print(f"\nRatio:{Ratio_transmitted_untransmitted}\n")
        mt2 = mt_trans_count.annotate_rows(
            variant_transmitted_singletons=hl.agg.count_where(mt_trans_count.transmitted_singletons_count == 1)
        )
        mt2.variant_transmitted_singletons.summarize()

        mt3 = mt_untrans_count.annotate_rows(
            variant_untransmitted_singletons=hl.agg.count_where(mt_untrans_count.untransmitted_singletons_count == 1)
        )
        mt3.variant_untransmitted_singletons.summarize()

        variant_transmitted_annot = mt2.rows()[ht.key].variant_transmitted_singletons
        variant_untransmitted_annot = mt3.rows()[ht.key].variant_untransmitted_singletons

    else:
        print("=== No pedifree data found: skipping transmitted/untransmitted calculations")
        variant_transmitted_annot = 0
        variant_untransmitted_annot = 0

    ht = ht.annotate(variant_transmitted_singletons=variant_transmitted_annot)
    ht = ht.annotate(variant_untransmitted_singletons=variant_untransmitted_annot)
    return ht


def transmitted_singleton_annotation(
    family_annot_htfile: str, trio_mtfile: str, trio_filtered_mtfile: str, trans_sing_htfile: str, mock: bool = False
) -> None:
    """
    Annotate MT with transmited singletons and create final variant QC HT for ranking
    :param str family_annot_htfile: Family annotation hail table file
    :param str trio_mtfile: Trio annotation hail matrixtable file
    :param str trio_filtered_mtfile: Trio filtered hail matrixtable file
    :param str trans_sing_htfile: Variant QC hail table file with transmitted singleton annotation
    """
    print("=== Annotating with transmitted singletons")
    ht = hl.read_table(family_annot_htfile)

    if not mock:
        # Open and filter the trios table
        mt_trios = hl.read_matrix_table(trio_mtfile)
        mt_trios = mt_trios.annotate_rows(consequence=ht[mt_trios.row_key].consequence)

        mt_filtered = mt_trios.filter_rows(
            (mt_trios.varqc_trios.AC[1] <= 2) & (mt_trios.consequence == "synonymous_variant")
        )
        mt_filtered.write(trio_filtered_mtfile, overwrite=True)

    ht = count_trans_untransmitted_singletons(mt_filtered if not mock else None, ht)
    ht.write(trans_sing_htfile, overwrite=True)


# def run_tdt(mtfile: str, trans_sing_htfile: str, pedfile: str, tdt_htfile: str):
#     '''
#     Run transmission disequilibrium test to get counts of number transmitted and number unstransmitted for each
#     variant and annotate the RF output file with this
#     :param str mtfile: Mtfile used to create random forest input table
#     :param str trans_sing_htfile: Variant QC hail table file with transmitted singleton annotation
#     :param str pedfile: Pedfile
#     :param str tdt_htfile: Htfile with transmitted/untransmitted counts
#     '''
#     print("Annotating with transmitted/unstransmitted counts")
#     mt = hl.read_matrix_table(mtfile)
#     pedigree = hl.Pedigree.read(pedfile)
#     tdt_ht = hl.transmission_disequilibrium_test(mt, pedigree)
#     ht = hl.read_table(trans_sing_htfile)
#     ht=ht.annotate(n_transmitted=tdt_ht[ht.key].t)
#     ht=ht.annotate(n_untransmitted=tdt_ht[ht.key].u)
#     ht.write(tdt_htfile, overwrite = True)


def annotate_gnomad(tdt_htfile: str, gnomad_htfile: str, final_htfile: str) -> None:
    """
    Annotate with gnomad allele frequencies
    :param str tdt_htfile: Htfile with transmitted/untransmitted counts
    :param str gnomad_htfile: Gnomad annotation hail table file
    :param str final_htfile: Final RF htfile for ranking and binning
    """
    print("=== Annotating with gnomad AF")
    ht = hl.read_table(tdt_htfile)
    gnomad_ht = hl.read_table(gnomad_htfile)
    ht = ht.annotate(gnomad_af=gnomad_ht[ht.key].freq[0].AF)
    ht.write(final_htfile, overwrite=True)


def get_options() -> Any:
    """
    Get options from the command line
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--syn-cqs", help="Annotate with synonymous CQs", action="store_true")
    parser.add_argument("--family-stats", help="Annotate with family stats and DNMs", action="store_true")
    parser.add_argument("--transmitted", help="Annotate with transmitted singletons", action="store_true")
    parser.add_argument("--gnomad-af", help="Annotate with gnomad AF", action="store_true")
    parser.add_argument("--all", help="run all steps", action="store_true")
    args = parser.parse_args()

    return args


def main() -> None:
    # set up
    args = get_options()
    inputs = parse_config()

    data_root: str = inputs["data_root"]
    dataset_name: str = inputs["dataset_name"]
    mtdir: str = os.path.join(data_root, inputs["matrixtables_lustre_dir"])
    resourcedir = os.path.join(data_root, inputs["resource_dir"])
    annot_dir: str = os.path.join(data_root, inputs["annotation_lustre_dir"])
    rf_dir = os.path.join(data_root, inputs["var_qc_rf_dir"])
    run_hash = inputs["runhash"]
    ped_file_name = inputs["pedfile_name"]

    # initialise hail
    sc = hail_utils.init_hl(inputs["tmp_dir"])

    htfile = os.path.join(rf_dir, run_hash, "rf_result.ht")
    synonymous_file = os.path.join(annot_dir, "synonymous_variants.txt")
    ht_cq_file = os.path.join(rf_dir, run_hash, "rf_result_with_synonymous.ht")
    if args.syn_cqs or args.all:
        print("=== Annotating with synonymous CQs")
        add_cq_annotation("file://" + htfile, "file://" + synonymous_file, "file://" + ht_cq_file)

    dnm_htfile = os.path.join(mtdir, "denovo_table.ht")
    fam_stats_htfile = os.path.join(mtdir, "family_stats.ht")
    trio_stats_htfile = os.path.join(mtdir, "trio_stats.ht")
    family_annot_htfile = os.path.join(rf_dir, run_hash, "rf_result_denovo_family_stats.ht")
    if args.family_stats or args.all:
        print("=== Annotate with family stats and DNMs")
        dnm_and_family_annotation(
            "file://" + ht_cq_file,
            "file://" + dnm_htfile,
            "file://" + fam_stats_htfile,
            "file://" + trio_stats_htfile,
            "file://" + family_annot_htfile,
            mock=(ped_file_name == ""),
        )

    trio_mtfile = os.path.join(mtdir, "trios.mt")
    trio_filtered_mtfile = os.path.join(mtdir, "trios_filtered_.mt")
    trans_sing_htfile = os.path.join(rf_dir, run_hash, "rf_result_trans_sing.ht")
    if args.transmitted or args.all:
        print("=== Annotate with transmitted singletons")
        transmitted_singleton_annotation(
            "file://" + family_annot_htfile,
            "file://" + trio_mtfile,
            "file://" + trio_filtered_mtfile,
            "file://" + trans_sing_htfile,
            mock=(ped_file_name == ""),
        )

    final_htfile = os.path.join(rf_dir, run_hash, "rf_result_final_for_ranking.ht")
    gnomad_htfile = os.path.join(resourcedir, "gnomad.exomes.r2.1.1.sites.liftover_grch38.ht")
    if args.gnomad_af or args.all:
        print("=== Annotate with gnomad AF")
        annotate_gnomad("file://" + trans_sing_htfile, "file://" + gnomad_htfile, "file://" + final_htfile)

    hail_utils.stop_hl(sc)


if __name__ == "__main__":
    main()
