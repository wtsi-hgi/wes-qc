# create hail table for random forest
import os

import hail as hl
import utils.constants as constants
from utils.utils import parse_config
from gnomad.variant_qc.random_forest import median_impute_features

from wes_qc import hail_utils


def create_rf_ht(
    mtfile: str,
    truthset_file: str,
    trio_stats_file: str,
    allele_data_file: str,
    allele_counts_file: str,
    inbreeding_file: str,
    htfile_rf_all_cols: str,
    htfile_rf_var_type_all_cols: str,
) -> None:
    """
    Load input mt and training data to create an input for random forest
    param str mtfile: Input matrixtable file
    param str truthset_file: Truthset hail table file
    param str trio_stats_file: Trio stats hail table file
    param str allele_data_file: Allele data hail table file
    param str allele_counts_file: Allele counts hail table file
    param str inbreeding_file: Inbreeding hail table file
    param str htfile_rf_all_cols: Output file for RF hail table
    param str htfile_rf_var_type_all_cols: Output file for RF hail table by variant type
    """
    n_partitions = 200
    mt = hl.read_matrix_table(mtfile)
    truth_data_ht = hl.read_table(truthset_file)
    allele_data_ht = hl.read_table(allele_data_file)
    allele_counts_ht = hl.read_table(allele_counts_file)
    inbreeding_ht = hl.read_table(inbreeding_file)

    allele_counts_ht = allele_counts_ht.select(*["ac_qc_samples_raw", "ac_qc_samples_adj"])

    # mt = mt.key_rows_by('locus').distinct_by_row().key_rows_by('locus', 'alleles')
    mt = mt.select_entries(GT=hl.unphased_diploid_gt_index_call(mt.GT.n_alt_alleles()))
    mt = mt.annotate_rows(InbreedingCoeff=hl.or_missing(~hl.is_nan(mt.info.InbreedingCoeff), mt.info.InbreedingCoeff))

    ht = mt.rows()
    ht = ht.transmute(**ht.info)
    ht = ht.select("MQ", "InbreedingCoeff", "a_index", "was_split", "meanHetAB", *constants.INFO_FEATURES)
    ht = ht.annotate(
        **inbreeding_ht[ht.key],
        **truth_data_ht[ht.key],
        **allele_data_ht[ht.key].allele_data,
        **allele_counts_ht[ht.key],
    )

    group = "raw"
    if trio_stats_file != "":  # Annotating the whole table with the trios information
        trio_stats_table = hl.read_table(trio_stats_file)
        trio_stats_ht = trio_stats_table.select(f"n_transmitted_{group}", f"ac_children_{group}")
        ht = ht.annotate(**trio_stats_ht[ht.key])
    else:  # Making fake annotations for non-present trios
        fake_trio_dict = {f"n_transmitted_{group}": 0, f"ac_children_{group}": 0}
        ht = ht.annotate(**fake_trio_dict)

    # annotate with C>A or not
    ht = ht.annotate(
        is_CA=((ht.alleles[0] == "C") & (ht.alleles[1] == "A")) | ((ht.alleles[0] == "G") & (ht.alleles[1] == "T"))
    )
    # annotate with all other possible SNPs
    # ht = ht.annotate(is_AC=((ht.alleles[0] == "A") & (ht.alleles[1] == "C")) | ((ht.alleles[0] == "T") & (ht.alleles[1] == "G")))
    # ht = ht.annotate(is_AG=((ht.alleles[0] == "A") & (ht.alleles[1] == "G")) | ((ht.alleles[0] == "T") & (ht.alleles[1] == "C")))
    # ht = ht.annotate(is_AT=((ht.alleles[0] == "A") & (ht.alleles[1] == "T")) | ((ht.alleles[0] == "T") & (ht.alleles[1] == "A")))
    # ht = ht.annotate(is_CG=((ht.alleles[0] == "C") & (ht.alleles[1] == "G")) | ((ht.alleles[0] == "G") & (ht.alleles[1] == "C")))
    # ht = ht.annotate(is_CT=((ht.alleles[0] == "C") & (ht.alleles[1] == "T")) | ((ht.alleles[0] == "G") & (ht.alleles[1] == "A")))

    ht = ht.annotate(fail_hard_filters=(ht.QD < 2) | (ht.FS > 60) | (ht.MQ < 30))
    ht = ht.annotate(ac_raw=ht.ac_qc_samples_raw)

    ht = ht.annotate(transmitted_singleton=(ht[f"n_transmitted_{group}"] == 1) & (ht[f"ac_qc_samples_{group}"] == 2))

    ht = ht.select(
        "a_index",
        # "was_split",
        *constants.FEATURES,
        *constants.TRUTH_DATA,
        **{
            "transmitted_singleton": (ht[f"n_transmitted_{group}"] == 1) & (ht[f"ac_qc_samples_{group}"] == 2),
            "fail_hard_filters": (ht.QD < 2) | (ht.FS > 60) | (ht.MQ < 30),
        },
        ac_raw=ht.ac_qc_samples_raw,
    )

    ht = ht.repartition(n_partitions, shuffle=False)
    ht = ht.checkpoint(htfile_rf_all_cols, overwrite=True)
    ht = median_impute_features(ht, strata={"variant_type": ht.variant_type})
    ht = ht.checkpoint(htfile_rf_var_type_all_cols, overwrite=True)


def main() -> None:
    # set up
    inputs = parse_config()

    data_root: str = inputs["data_root"]
    dataset_name: str = inputs["dataset_name"]
    mtdir: str = os.path.join(data_root, inputs["matrixtables_lustre_dir"])
    resourcedir = os.path.join(data_root, inputs["resource_dir"])
    annot_dir: str = os.path.join(data_root, inputs["annotation_lustre_dir"])

    training_sets_dir = os.path.join(data_root, inputs["training_set_dir"])
    ped_file_name = inputs["pedfile_name"]

    # initialise hail
    sc = hail_utils.init_hl(inputs["tmp_dir"])

    truthset_file = os.path.join(resourcedir, "truthset_table.ht")
    trio_stats_file = os.path.join(mtdir, "trio_stats.ht")
    allele_data_file = os.path.join(mtdir, "allele_data.ht")
    allele_counts_file = os.path.join(mtdir, "qc_ac.ht")
    inbreeding_file = os.path.join(mtdir, "inbreeding.ht")

    # mtfile = mtdir + "mt_pops_QC_filters_sequencing_location_and_superpop_sanger_only_after_sample_qc.mt"
    mtfile = os.path.join(mtdir, "mt_varqc_splitmulti.mt")
    htfile_rf_all_cols = os.path.join(mtdir, "ht_for_RF_all_cols.ht")
    htfile_rf_var_type_all_cols = os.path.join(mtdir, "ht_for_RF_by_variant_type_all_cols.ht")

    create_rf_ht(
        "file://" + mtfile,
        "file://" + truthset_file,
        ("file://" + trio_stats_file) if ped_file_name != "" else "",
        "file://" + allele_data_file,
        "file://" + allele_counts_file,
        "file://" + inbreeding_file,
        "file://" + htfile_rf_all_cols,
        "file://" + htfile_rf_var_type_all_cols,
    )

    hail_utils.stop_hl(sc)


if __name__ == "__main__":
    main()
