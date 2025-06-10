# create hail table for random forest
from typing import Optional

import hail as hl
import utils.constants as constants
from utils.utils import parse_config, path_spark
from gnomad.variant_qc.random_forest import median_impute_features

from wes_qc import hail_utils


def prepare_matrix_table(mt: hl.MatrixTable) -> hl.MatrixTable:
    """Prepare the matrix table by setting GT and InbreedingCoeff."""
    mt = mt.select_entries(GT=hl.unphased_diploid_gt_index_call(mt.GT.n_alt_alleles()))
    mt = mt.annotate_rows(InbreedingCoeff=hl.or_missing(~hl.is_nan(mt.info.InbreedingCoeff), mt.info.InbreedingCoeff))
    return mt


def create_initial_table(mt: hl.MatrixTable) -> hl.Table:
    """Create initial Hail table with basic features."""
    ht = mt.rows()
    ht = ht.transmute(**ht.info)
    # TODO: add workaround when AS_QD and related features are not accessible in the initial dataset
    ht = ht.select("MQ", "InbreedingCoeff", "a_index", "was_split", "meanHetAB", *constants.INFO_FEATURES)
    return ht


def annotate_with_external_data(
    ht: hl.Table, inbreeding_ht: hl.Table, truth_data_ht: hl.Table, allele_data_ht: hl.Table, allele_counts_ht: hl.Table
) -> hl.Table:
    """Annotate table with external data."""
    return ht.annotate(
        **inbreeding_ht[ht.key],
        **truth_data_ht[ht.key],
        **allele_data_ht[ht.key].allele_data,
        **allele_counts_ht[ht.key],
    )


def add_trio_stats(
    ht: hl.Table,
    trio_stats_table: Optional[hl.Table],
    group: str = "raw",  # TODO: do we need this parameter
) -> hl.Table:
    """Add trio statistics to the table."""
    if trio_stats_table is not None:
        print("--- Annotation with trio stats ---")
        trio_stats_ht = trio_stats_table.select(f"n_transmitted_{group}", f"ac_children_{group}")
        return ht.annotate(**trio_stats_ht[ht.key])

    print("--- No trio stats provided - skipping trio annotation ---")
    fake_trio_dict = {f"n_transmitted_{group}": 0, f"ac_children_{group}": 0}
    return ht.annotate(**fake_trio_dict)


def add_filters_and_singletons(
    ht: hl.Table,
    fail_hard_filters_QD_less_than: float,
    fail_hard_filters_FS_greater_than: float,
    fail_hard_filters_MQ_less_than: float,
    group: str = "raw",
) -> hl.Table:
    """Add hard filters and singleton annotations."""
    ht = ht.annotate(
        fail_hard_filters=(ht.QD < fail_hard_filters_QD_less_than)
        | (ht.FS > fail_hard_filters_FS_greater_than)
        | (ht.MQ < fail_hard_filters_MQ_less_than)
    )
    ht = ht.annotate(ac_raw=ht.ac_qc_samples_raw)
    return ht.annotate(transmitted_singleton=(ht[f"n_transmitted_{group}"] == 1) & (ht[f"ac_qc_samples_{group}"] == 2))


def select_final_features(
    ht: hl.Table,
    fail_hard_filters_QD_less_than: float,
    fail_hard_filters_FS_greater_than: float,
    fail_hard_filters_MQ_less_than: float,
    group: str = "raw",
) -> hl.Table:
    """Select final features for the random forest."""
    # TODO: are these thresholds separate from the ones above?

    return ht.select(
        "a_index",
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
    **kwargs,
) -> hl.Table:
    """
    Load input mt and training data to create an input for random forest
    """
    n_partitions = 200
    group = "raw"
    # Prepare input data
    allele_counts_ht = allele_counts_ht.select(*["ac_qc_samples_raw", "ac_qc_samples_adj"])

    mt = prepare_matrix_table(mt)

    # Create and annotate table
    ht = create_initial_table(mt)
    ht = annotate_with_external_data(ht, inbreeding_ht, truth_data_ht, allele_data_ht, allele_counts_ht)
    ht = add_trio_stats(ht, trio_stats_table, group)

    ht = add_filters_and_singletons(
        ht, fail_hard_filters_QD_less_than, fail_hard_filters_FS_greater_than, fail_hard_filters_MQ_less_than, group
    )
    ht = select_final_features(
        ht, fail_hard_filters_QD_less_than, fail_hard_filters_FS_greater_than, fail_hard_filters_MQ_less_than, group
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
    print(
        "pedfile: ",
        pedfile,
    )
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
