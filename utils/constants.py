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
    "InbreedingCoeff",
    "variant_type",
    "allele_type",
    "n_alt_alleles",
    "was_mixed",
    "has_star",
    "MQ",
    "QD",
    "MQRankSum",
    "SOR",
    "ReadPosRankSum",
]

TRUTH_DATA = ["hapmap", "omni", "mills", "kgp_phase1_hc"]