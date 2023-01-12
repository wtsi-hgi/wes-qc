#constants for random forest

INFO_FEATURES = [
    "AQ"
] 

FEATURES = [
    "variant_type",
    "allele_type",
    "n_alt_alleles",
    "was_mixed",
    "was_split",
    "has_star",
    "AQ_allele",
    "is_CA",
    "meanHetAB"
]

TRUTH_DATA = ["hapmap", "omni", "mills", "kgp_phase1_hc"]

LABEL_COL = "rf_label"
TRAIN_COL = "rf_train"
PREDICTION_COL = "rf_prediction"