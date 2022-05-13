#create hail table for random forest
import hail as hl
import pyspark
import wes_qc.utils.constants as constants
from wes_qc.utils.utils import parse_config
from gnomad.variant_qc.random_forest import median_impute_features


# INFO_FEATURES = [
#     "AS_QD",
#     "AS_ReadPosRankSum",
#     "AS_MQRankSum",
#     "AS_SOR",
#     "QD",
#     "MQRankSum",
#     "SOR",
#     "ReadPosRankSum",
#     "FS",
#     "DP"
# ] 

# FEATURES = [
#     "InbreedingCoeff",
#     "variant_type",
#     "allele_type",
#     "n_alt_alleles",
#     "was_mixed",
#     "has_star",
#     "MQ",
#     "QD",
#     "MQRankSum",
#     "SOR",
#     "ReadPosRankSum",
# ]

# TRUTH_DATA = ["hapmap", "omni", "mills", "kgp_phase1_hc"]


def create_rf_ht(mtfile: str, truthset_file: str, trio_stats_file: str, allele_data_file: str, allele_counts_file: str, inbreeding_file: str, htfile_rf_all_cols: str, htfile_rf_var_type_all_cols: str):
    '''
    Load input mt and training data to create an input for random forest
    param str mtfile: Input matrixtable file
    param str truthset_file: Truthset hail table file
    param str trio_stats_file: Trio stats hail table file
    param str allele_data_file: Allele data hail table file
    param str allele_counts_file: Allele counts hail table file
    param str inbreeding_file: Inbreeding hail table file
    param str htfile_rf_all_cols: Output file for RF hail table 
    param str htfile_rf_var_type_all_cols: Output file for RF hail table by variant type
    '''
    n_partitions = 200
    mt = hl.read_matrix_table(mtfile)
    truth_data_ht = hl.read_table(truthset_file)
    trio_stats_table = hl.read_table(trio_stats_file)
    allele_data_ht = hl.read_table(allele_data_file)
    allele_counts_ht = hl.read_table(allele_counts_file)
    inbreeding_ht = hl.read_table(inbreeding_file)

    allele_counts_ht = allele_counts_ht.select(*['ac_qc_samples_raw', 'ac_qc_samples_adj'])
    group = "raw"
    trio_stats_ht = trio_stats_table.select(
        f"n_transmitted_{group}", f"ac_children_{group}"
    )

    mt = mt.key_rows_by('locus').distinct_by_row().key_rows_by('locus', 'alleles')
    mt = mt.select_entries(GT=hl.unphased_diploid_gt_index_call(mt.GT.n_alt_alleles()))
    mt = mt.annotate_rows(InbreedingCoeff=hl.or_missing(~hl.is_nan(mt.info.InbreedingCoeff), mt.info.InbreedingCoeff))

    ht = mt.rows()
    ht = ht.transmute(**ht.info)
    ht = ht.select( "MQ", "InbreedingCoeff", *constants.INFO_FEATURES)
    ht = ht.annotate(
        **inbreeding_ht[ht.key],
        **trio_stats_ht[ht.key],
        **truth_data_ht[ht.key],
        **allele_data_ht[ht.key].allele_data,
        **allele_counts_ht[ht.key],
    )

    ht = ht.annotate(fail_hard_filters=(ht.QD < 2)
                     | (ht.FS > 60) | (ht.MQ < 30))
    ht = ht.annotate(ac_raw=ht.ac_qc_samples_raw)
    ht = ht.annotate(transmitted_singleton=(
        ht[f"n_transmitted_{group}"] == 1) & (ht[f"ac_qc_samples_{group}"] == 2))

    ht = ht.select(
        "a_index",
        "was_split",
        *constants.FEATURES,
        *constants.TRUTH_DATA,
        **{
            "transmitted_singleton": (ht[f"n_transmitted_{group}"] == 1)
            & (ht[f"ac_qc_samples_{group}"] == 2),
            "fail_hard_filters": (ht.QD < 2) | (ht.FS > 60) | (ht.MQ < 30),
        },
        ac_raw=ht.ac_qc_samples_raw
    )

    ht = ht.repartition(n_partitions, shuffle=False)
    ht = ht.checkpoint(htfile_rf_all_cols, overwrite=True)
    ht = median_impute_features(ht, strata={"variant_type": ht.variant_type})
    ht = ht.checkpoint(htfile_rf_var_type_all_cols, overwrite=True)


def main():
    #set up
    inputs = parse_config()
    mtdir = inputs['matrixtables_lustre_dir']
    resourcedir = inputs['resource_dir']

    # initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    truthset_file = resourcedir + "truthset_table.ht"
    trio_stats_file = mtdir + "trio_stats.ht"
    allele_data_file = mtdir + "allele_data.ht"
    allele_counts_file = mtdir + "qc_ac.ht"
    inbreeding_file = mtdir + "inbreeding.ht"

    mtfile = mtdir + "mt_pops_QC_filters_sequencing_location_and_superpop_sanger_only_after_sample_qc.mt"
    htfile_rf_all_cols = mtdir + "ht_for_RF_all_cols.ht"
    htfile_rf_var_type_all_cols = mtdir + "ht_for_RF_by_variant_type_all_cols.ht"

    create_rf_ht(mtfile, truthset_file, trio_stats_file, allele_data_file, allele_counts_file, inbreeding_file, htfile_rf_all_cols, htfile_rf_var_type_all_cols)


if __name__ == '__main__':
    main()