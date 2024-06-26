#constants for random forest

INFO_FEATURES = [
    "AS_QD",
    "AS_ReadPosRankSum",
    "AS_MQRankSum",
    "AS_SOR",
    "QD",
    "MQRankSum",
    "SOR",
    "ReadPosRankSum",
    "FS",
    "DP"
] 

FEATURES = [
    "variant_type",
    "allele_type",
    "n_alt_alleles",
    "was_mixed",
    "was_split",
    "has_star",
    "MQ",
    "QD",
    "MQRankSum",
    "SOR",
    "ReadPosRankSum",
    "meanHetAB"
]

TRUTH_DATA = ["hapmap", "omni", "mills", "kgp_phase1_hc"]

LABEL_COL = "rf_label"
TRAIN_COL = "rf_train"
PREDICTION_COL = "rf_prediction"