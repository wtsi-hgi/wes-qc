---
tags:
  - HowTo
  - Compute
  - OpenStack
  - Hail
  - Pipelines
---

# Getting Started With WES QC Using Hail

!!! warning
    This documentation is under development and may be incomplete.

This guide covers WES QC using Hail. It is important to note that every dataset is different and that for best results it is not advisible to view this guide as a recipie for QC. Each dataset will require careful tailoring and evaluation of the QC for best results.

## Before you start

In order to run through this guide you will need an openstack cluster with Hail and Spark installed. It is recommended that you use osdataproc to create this.

This guide also requires a WES dataset joint called with GATK using the cromwell pipeline. If starting with a Hail matrixtable then start at stage 2, sample QC.

## Set up

Each dataset should have its own branch of the WES QC code. Clone the main branch using:
`git clone https://gitlab.internal.sanger.ac.uk/hgi-projects/birth_cohort_wes/wes_qc.git`

or

`git clone git@gitlab.internal.sanger.ac.uk:hgi-projects/birth_cohort_wes/wes_qc.git`

Create a new branch and switch to your branch:

`cd wes_qc`

`git checkout -b my_project`

Edit *config/inputs.yaml* to include the correct paths for your datsets and working directories. You will need to edit al directory paths in this file and commit/push your changes to git.

Start a new tmux session, and edit the PYTHONPATH to include the directory you originally cloned the git repo into.

`tmux`

`export PYTHONPATH=/home/ubuntu/:$PYTHONPATH`

### 1. Load data fronm VCFs into Hail

Copy the VCFs from the cromwell directories (cromwell_base_dir and regions in the config) into an area accessible via your hail clusater (target dir in the config). This step can be done on the farm

`1-import_data/1-copy_gatk_to_secure_lustre.py`

Subsequent steps are run on the Hail cluster using park-submit.

Load VCFs into Hail and save as a Hail MatrixTable

`PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 1-import_data/2-import_gatk_vcfs_to_hail.py`

### 2. Sample QC

Apply hard filters and annotate with imputed sex

`PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 2-sample_qc/1-hard_filters_sex_annotation.py`

Prune by linkage disequilibrium, and identify samples from related individuals with PCrelate

`PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 2-sample_qc/2-prune_related_samples.py`

The following five steps deal with prediction of superpopulation by PCA and projection on to 1000 genomes data of known population. If a matrixtable of 1000 genomes data is available the first step can be omitted.

Create MatrixTable from 1kg WES VCFs

`PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 2-sample_qc/3-population_pca_prediction.py --kg_to_mt`

Merge 1kg MatrixTable with WES MatrixTable

`PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 2-sample_qc/3-population_pca_prediction.py --merge`

Annotate with known populations for 1kg samples and filter to remove low quality variants, variants in long range linkage disequilibrium regions and palindromic SNPs

`PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 2-sample_qc/3-population_pca_prediction.py --filter`

Run PCA

`PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 2-sample_qc/3-population_pca_prediction.py --pca`

Run population prediction

`PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 2-sample_qc/3-population_pca_prediction.py --assign_pops`

Use a Jupyter Notebook on the Hail cluster to plot PC1 vs PC2 and PC1 vs PC3 for known populations (1000 genomes only) and predicted populations. An example Jupyter Notebook is in *jupyter_notebooks/sample_qc/population_inference_plots.ipynb*

Now that we have the predicted populations that each sample belongs to we can run sample QC stratified by population and identify outliers within each population for each metric tested: Number of SNPs, number of deletions, number of insertions, insertion/deletion ratip, heterozygosity rate, heterozygous/homozygous ratio, number of transitions, number of transversions, transition/transversion ratio.

`PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 2-sample_qc/4-find_population_outliers.py`

Filter data to exclude samples which fail QC (note step 5 in the git repository is for non-stratified sample QC and is not included in this tutorial). The final step in sample QC is filtering the data to remove samples which are identifed as failing in the previous script. At this stage samples failing on FREEMIX score and on identity checks are also removed. This samples should be in files in the annotations directory: verify_bam_id_result_concat.selfSM lists sample ID and FREEMIX score and sanger_samples_excluded_after_gtcheck.txt lists samples failing identity checks. If no samples fail identify checks the latter file could be empty.

`PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 2-sample_qc/6-filter_fail_sample_qc.py`

### 3. Variant QC

Variant QC requires a ped file detailing any trios. This is an unheadered tab delinmtied file containing the following columns:
- Family ID
- Proband ID
- Paternal ID
- Maternal ID
- Proband sex (1 or 2)
- Proband affected status (0 or 1)

The first step of variant QC produces a truth set of known true positive variants and generates annotation.

`PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 3-variant_qc/1-generate_truth_sets.py --all`

Next an input table is generated to run the random forest on

`PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 3-variant_qc/2-create_rf_ht.py`

Next train the random forest model, note the use of --master local[*] here which is needed to ensure that the worker nodes pick up the correct python modules.

`PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit --master local[*] 3-variant_qc/3-train_rf.py`

A run hash ID will be printed to STDOUT, this is needed for future steps, so note it down. It is an 8 character string consisting of letters and numbers and is represented in the following commands as run_hash.

Now apply the random forest to the entire dataset.

`PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit --master local[*] 3-variant_qc/4-apply_rf.py --runhash run_hash`

Annotate the random forest output with metrics including synonymous variants, family annotation, transmitted/untransmitted singletons and gnomAD allele frequency. Synonymous variants are required in a file generated from VEP annotation and in the following format:

```
chr10   100202145   rs200461553 T   G   synonymous_variant  
chr10   100204510   rs2862988   C   T   synonymous_variant  
chr10   100204528   rs374991603 G   A   synonymous_variant  
chr10   100204555   rs17880383  G   A   synonymous_variant  
```

`PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 3-variant_qc/5-annotate_ht_after_rf.py  --runhash run_hash`

Add ranks to variants based on random forest score, and bin the variants based on this.

`PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 3-variant_qc/6-rank_and_bin.py --runhash run_hash`

Create plots of the binned random forest output to use in the selection of thresholds. Separate thresholds are used for SNPs and indels.

`PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 3-variant_qc/7-plot_rf_output.py --runhash run_hash`

After examnining the plots and selecting suitable thresholds for SNPs and indels you can calculate the number of true positive and false positive variants remaining at your chosen thresholds and at the bins surrounding it using the following (where snp_bin and indel_bin are the thresholds selected for SNPs and indels respectively).

`PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 3-variant_qc/8-select_thresholds.py --runhash run_hash --snv snp_bin --indel indel_bin`

Filter the variants in the Hail MatrixTable based on the selected threshold for SNPs and indels

`PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 3-variant_qc/9-filter_mt_after_variant_qc.py --runhash run_hash --snv snp_bin --indel indel_bin`

### 4. Genotype QC

Genotype QC may be performed using a range of filters defined in *config/inputs.yaml* or it may be performed using a single set of filters for DP, GQ, VAF and bins. Here we apply a range of filters (relaxed, medium and stringent) which we have defined in the *config/inputs.yaml* file

`PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 4-genotype_qc/1a-apply_range_of_hard_filters.py --runhash run_hash`

In order to evaluate these filters, variant counts per consequence per sample and transmitted/untransmitted ratio of synonymous singletons (if trios are present in the data) are calculated as follows. VEP annotation is required for this step in the following format:

```
chr10   100199947   rs367984062 A   C   intron_variant  
chr10   100199976   rs774723210 G   A   missense_variant  
chr10   100200004   .   C   A   missense_variant  
chr10   100200012   rs144642900 C   T   missense_variant  
chr1    100200019   .  C    A   stop_gained&splice_region_variant  
```

`PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 4-genotype_qc/2-counts_per_sample.py`

Export the filtered variants to VCF

`PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit 4-genotype_qc/3-export_vcfs.py`


