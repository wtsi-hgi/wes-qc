#!/bin/bash

# Set up
export PYTHONPATH=/lustre/scratch126/dh24_test/wes-qc
export WES_CONFIG=/lustre/scratch126/dh24_test/wes-qc/config/inputs-mk43-a.yaml


# Step 1
PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 1-import_data/1-import_gatk_vcfs_to_hail.py

# Step 2
PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 2-sample_qc/1-hard_filters_sex_annotation.py
PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 2-sample_qc/2-prune_related_samples.py
PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 2-sample_qc/3-population_pca_prediction.py --kg_to_mt
PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 2-sample_qc/3-population_pca_prediction.py --run
PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 2-sample_qc/4-find_population_outliers.py
PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 2-sample_qc/5-filter_fail_sample_qc.py

# Step 3
PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 3-variant_qc/variant_qc_non_trios/1-generate_truth_sets_non_trios.py --all
PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 3-variant_qc/variant_qc_non_trios/2-create_rf_ht_non_trios.py
PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit --master local[*] 3-variant_qc/3-train_rf.py --manual_runhash controllable_test_run_id
PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit --master local[*] 3-variant_qc/4-apply_rf.py --runhash controllable_test_run_id
PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 3-variant_qc/variant_qc_non_trios/5-annotate_ht_after_rf_no_trios.py --runhash controllable_test_run_id
PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 3-variant_qc/variant_qc_non_trios/6-rank_and_bin_no_trios.py --runhash controllable_test_run_id
PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 3-variant_qc/variant_qc_non_trios/7-plot_rf_output_no_trios.py --runhash controllable_test_run_id
PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 3-variant_qc/8-select_thresholds.py --runhash controllable_test_run_id --snv 92 --indel 68
PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 3-variant_qc/9-filter_mt_after_variant_qc.py --runhash controllable_test_run_id --snv 84 --indel 60

# Step 4
PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 4-genotype_qc/1-apply_hard_filters.py --dp 5 --gq 10 --ab 0.2
PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 4-genotype_qc/1a-apply_range_of_hard_filters.py --runhash controllable_test_run_id
PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 4-genotype_qc/2-counts_per_sample.py
PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 4-genotype_qc/3-export_vcfs.py
PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 4-genotype_qc/3a-export_vcfs_range_of_hard_filters.py --runhash controllable_test_run_id
PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 4-genotype_qc/3b-export_vcfs_stingent_filters.py --runhash controllable_test_run_id

