# create hail table for random forest
import hail as hl
import utils.constants as constants
from utils.utils import parse_config, path_spark
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
    config: dict,
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
    conf = config["step3"]["create_rf_ht"]

    n_partitions = 200
    mt = hl.read_matrix_table(path_spark(mtfile))
    truth_data_ht = hl.read_table(path_spark(truthset_file))
    trio_stats_table = hl.read_table(path_spark(trio_stats_file))
    allele_data_ht = hl.read_table(path_spark(allele_data_file))
    allele_counts_ht = hl.read_table(path_spark(allele_counts_file))
    inbreeding_ht = hl.read_table(path_spark(inbreeding_file))

    allele_counts_ht = allele_counts_ht.select(*["ac_qc_samples_raw", "ac_qc_samples_adj"])
    group = "raw"
    trio_stats_ht = trio_stats_table.select(f"n_transmitted_{group}", f"ac_children_{group}")

    # mt = mt.key_rows_by('locus').distinct_by_row().key_rows_by('locus', 'alleles')
    mt = mt.select_entries(GT=hl.unphased_diploid_gt_index_call(mt.GT.n_alt_alleles()))
    mt = mt.annotate_rows(InbreedingCoeff=hl.or_missing(~hl.is_nan(mt.info.InbreedingCoeff), mt.info.InbreedingCoeff))

    ht = mt.rows()
    ht = ht.transmute(**ht.info)
    ht = ht.select("MQ", "InbreedingCoeff", "a_index", "was_split", "meanHetAB", *constants.INFO_FEATURES)
    ht = ht.annotate(
        **inbreeding_ht[ht.key],
        **trio_stats_ht[ht.key],
        **truth_data_ht[ht.key],
        **allele_data_ht[ht.key].allele_data,
        **allele_counts_ht[ht.key],
    )

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

    ht = ht.annotate(
        fail_hard_filters=(ht.QD < conf["fail_hard_filters_QD_less_than"])
        | (ht.FS > conf["fail_hard_filters_FS_greater_than"])
        | (ht.MQ < conf["fail_hard_filters_MQ_less_than"])
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
            "fail_hard_filters": (ht.QD < conf["fail_hard_filters_QD_less_than"])
            | (ht.FS > conf["fail_hard_filters_FS_greater_than"])
            | (ht.MQ < conf["fail_hard_filters_MQ_less_than"]),
        },
        ac_raw=ht.ac_qc_samples_raw,
    )

    ht = ht.repartition(n_partitions, shuffle=False)
    ht = ht.checkpoint(path_spark(htfile_rf_all_cols), overwrite=True)
    ht = median_impute_features(ht, strata={"variant_type": ht.variant_type})
    ht.write(path_spark(htfile_rf_var_type_all_cols), overwrite=True)


def main():
    # set up
    config = parse_config()
    tmp_dir = config["general"]["tmp_dir"]

    # = STEP LOGIC = #
    hail_utils.init_hl(tmp_dir)

    truthset_file = config["step3"]["create_rf_ht"]["truthset_file"]
    trio_stats_file = config["step3"]["create_rf_ht"]["trio_stats_file"]
    allele_data_file = config["step3"]["create_rf_ht"]["allele_data_file"]
    allele_counts_file = config["step3"]["create_rf_ht"]["allele_counts_file"]
    inbreeding_file = config["step3"]["create_rf_ht"]["inbreeding_file"]

    # mtfile = mtdir + "mt_pops_QC_filters_sequencing_location_and_superpop_sanger_only_after_sample_qc.mt"
    mtfile = config["step3"]["create_rf_ht"]["mtfile"]
    htoutfile_rf_all_cols = config["step3"]["create_rf_ht"]["htoutfile_rf_all_cols"]
    htoutfile_rf_var_type_all_cols = config["step3"]["create_rf_ht"]["htoutfile_rf_var_type_all_cols"]

    create_rf_ht(
        mtfile,
        truthset_file,
        trio_stats_file,
        allele_data_file,
        allele_counts_file,
        inbreeding_file,
        htoutfile_rf_all_cols,
        htoutfile_rf_var_type_all_cols,
        config,
    )


if __name__ == "__main__":
    main()
