# create hail table for random forest
from typing import Optional

import hail as hl
import utils.constants as constants
from utils.utils import parse_config, path_spark
from gnomad.variant_qc.random_forest import median_impute_features

from wes_qc import hail_utils


def create_rf_ht(
    mt: hl.MatrixTable,
    truth_data_ht: hl.Table,
    trio_stats_table: Optional[hl.Table],
    allele_data_ht: hl.Table,
    allele_counts_ht: hl.Table,
    inbreeding_ht: hl.Table,
    fail_hard_filters_QD_less_than: float,
    fail_hard_filters_FS_greater_than: float,
    fail_hard_filters_MQ_less_than: float,
    annotate_CA=False,
    **kwargs,
) -> hl.Table:
    """
    Load input mt and training data to create an input for random forest
    """
    n_partitions = 200

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
    if trio_stats_table is not None:  # Annotating the whole table with the trios information
        print("--- Annotation with trio stats ---")
        trio_stats_ht = trio_stats_table.select(f"n_transmitted_{group}", f"ac_children_{group}")
        ht = ht.annotate(**trio_stats_ht[ht.key])
    else:
        print("--- No trio stats provided - skipping tiro annotation ---")
        fake_trio_dict = {f"n_transmitted_{group}": 0, f"ac_children_{group}": 0}
        ht = ht.annotate(**fake_trio_dict)

    # TODO: do we need this for regular datasets? If yes, generalize
    # Optionally annotate with C>A based on configuration
    if annotate_CA:
        ht = ht.annotate(
            is_CA=((ht.alleles[0] == "C") & (ht.alleles[1] == "A")) | ((ht.alleles[0] == "G") & (ht.alleles[1] == "T"))
        )
        # annotate with all other possible SNPs
        # ht = ht.annotate(is_AC=((ht.alleles[0] == "A") & (ht.alleles[1] == "C")) | ((ht.alleles[0] == "T") & (ht.alleles[1] == "G")))
        # ht = ht.annotate(is_AG=((ht.alleles[0] == "A") & (ht.alleles[1] == "G")) | ((ht.alleles[0] == "T") & (ht.alleles[1] == "C")))
        # ht = ht.annotate(is_AT=((ht.alleles[0] == "A") & (ht.alleles[1] == "T")) | ((ht.alleles[0] == "T") & (ht.alleles[1] == "A")))
        # ht = ht.annotate(is_CG=((ht.alleles[0] == "C") & (ht.alleles[1] == "G")) | ((ht.alleles[0] == "G") & (ht.alleles[1] == "C")))
        # ht = ht.annotate(is_CT=((ht.alleles[0] == "C") & (ht.alleles[1] == "T")) | ((ht.alleles[0] == "G") & (ht.alleles[1] == "A")))

    ht = ht.annotate(
        fail_hard_filters=(ht.QD < fail_hard_filters_QD_less_than)
        | (ht.FS > fail_hard_filters_FS_greater_than)
        | (ht.MQ < fail_hard_filters_MQ_less_than)
    )
    ht = ht.annotate(ac_raw=ht.ac_qc_samples_raw)
    ht = ht.annotate(transmitted_singleton=(ht[f"n_transmitted_{group}"] == 1) & (ht[f"ac_qc_samples_{group}"] == 2))

    # TODO: are these thresholds separate from the ones above?
    ht = ht.select(
        "a_index",
        # "was_split",
        *constants.FEATURES,
        *constants.TRUTH_DATA,
        **{
            "transmitted_singleton": (ht[f"n_transmitted_{group}"] == 1) & (ht[f"ac_qc_samples_{group}"] == 2),
            "fail_hard_filters": (ht.QD < fail_hard_filters_QD_less_than)
            | (ht.FS > fail_hard_filters_FS_greater_than)
            | (ht.MQ < fail_hard_filters_MQ_less_than),
        },
        ac_raw=ht.ac_qc_samples_raw,
    )

    ht = ht.repartition(n_partitions, shuffle=False)
    # ht = ht.checkpoint(path_spark(htfile_rf_all_cols), overwrite=True) # TODO: check do we need this matrixtable. Remove if not used
    # TODO: doublecheck is it methodologically correct to impute features used for model training
    ht = median_impute_features(ht, strata={"variant_type": ht.variant_type})

    return ht


def main():
    # set up
    config = parse_config()
    tmp_dir = config["general"]["tmp_dir"]

    # = STEP PARAMETERS = #
    pedfile = config["step3"]["pedfile"]

    # = STEP DEPENDENCIES = #
    truthset_file = config["step3"]["create_rf_ht"]["truthset_file"]
    mtfile = config["step3"]["create_rf_ht"]["mtfile"]
    trio_stats_file = config["step3"]["create_rf_ht"]["trio_stats_file"]
    allele_data_file = config["step3"]["create_rf_ht"]["allele_data_file"]
    allele_counts_file = config["step3"]["create_rf_ht"]["allele_counts_file"]
    inbreeding_file = config["step3"]["create_rf_ht"]["inbreeding_file"]

    # = STEP OUTPUTS = #
    htoutfile_rf_var_type_all_cols = config["step3"]["create_rf_ht"]["htoutfile_rf_var_type_all_cols"]

    # = STEP LOGIC = #
    hail_utils.init_hl(tmp_dir)

    mt = hl.read_matrix_table(path_spark(mtfile))
    truth_data_ht = hl.read_table(path_spark(truthset_file))
    trio_stats_table = hl.read_table(path_spark(trio_stats_file)) if pedfile is not None else None
    allele_data_ht = hl.read_table(path_spark(allele_data_file))
    allele_counts_ht = hl.read_table(path_spark(allele_counts_file))
    inbreeding_ht = hl.read_table(path_spark(inbreeding_file))

    ht = create_rf_ht(
        mt,
        truth_data_ht,
        trio_stats_table,
        allele_data_ht,
        allele_counts_ht,
        inbreeding_ht,
        **config["step3"]["create_rf_ht"],
    )
    ht.write(path_spark(htoutfile_rf_var_type_all_cols), overwrite=True)


if __name__ == "__main__":
    main()
