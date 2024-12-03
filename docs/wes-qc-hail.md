# Getting Started With WES QC Using Hail

> [!WARNING]
> This documentation is under development and may be incomplete.

This guide covers WES QC using Hail. It is important to note that every dataset is different and that for best results it is not advisable to view this guide as a recipe for QC. Each dataset will require careful tailoring and evaluation of the QC for best results.

## Before you start

In order to run through this guide you will need an OpenStack cluster with Hail and Spark installed. It is recommended that you use `osdataproc` to create this. Follow the [[hail-on-spark]] guide to create such a cluster.

This guide also requires a WES dataset joint called with GATK using the cromwell pipeline. If starting with a Hail matrixtable then start at Step [[#2. Sample QC]].

## Set up

Clone the main branch using:
```shell
git clone https://github.com/wtsi-hgi/wes-qc.git
cd wes_qc
```

Historically, each dataset has its own branch, but now the HGI team is trying to merge all functionality into the main branch
and vary only the config file.
This process is not finished yet, so if you have to make any dataset-specific operations, you can do it in your own branch.

Create a new config file for your dataset.
By default, all scripts will use the config fine named `inputs.yaml`.
You can make a symlink for it to keep the config name meaningful.

```shell
cd config
cp public-dataset.yaml my_project.yaml
ln -snf my_project.yaml inputs.yaml
cd ..
```

Edit `config/my_project.yaml` to include the correct paths for your datasets and working directories.
Most probably you'll need to change only `dataset_name` and `data_root` entries.
All other paths are specified as relatives, so you won't need top edit it.

Start a new tmux session, and edit the PYTHONPATH to include the directory you originally cloned the git repo into.

```shell
export PYTHONPATH=/path/to/wes-qc
```

>[!NOTE]
>To submit jobs to SPARK you can either:
>```shell
># OPTION 1
>export PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python
># OPTION 2
>source /home/ubuntu/venv/bin/activate
>
>spark-submit /path/to/hail_script.py
>```

### 1. Load data from VCFs into Hail

Load VCFs into Hail and save as a Hail MatrixTable

```shell
spark-submit 1-import_data/1-import_gatk_vcfs_to_hail.py
```

Create the 1000G population prediction resource set.
This resource set is required for the super-population prediction on the population PCA step.

```shell
spark-submit 1-import_data/1-import_gatk_vcfs_to_hail.py --all
```

You need to create this resource set only once.
Then you can reuse it with any data cohort.


### 2. Sample QC

1. Apply hard filters and annotate with imputed sex

```shell
spark-submit 2-sample_qc/1-hard_filters_sex_annotation.py
```

2. Prune by linkage disequilibrium, and identify samples from related individuals with PCrelate

```shell
spark-submit 2-sample_qc/2-prune_related_samples.py
```
This step doesn't affect any other steps,
because, on the step 2/3 we're using the clustering approach with PCA scores projection,
However, this information can be useful to compare it with pedigree data and identify mislabeled samples.

3. Population prediction with PCA

Merge 1kg MatrixTable with WES MatrixTable and make LD pruning.

```shell
spark-submit 2-sample_qc/3-population_pca_prediction.py --merge-and-ldprune
```

Run PCA.

```shell
spark-submit 2-sample_qc/3-population_pca_prediction.py --pca
```

Plot 1KG PCA. On this step, all dataset samples should be labelled as `N/A`.

```shell
spark-submit 2-sample_qc/3-population_pca_prediction.py --pca-plot
```

Run population prediction.

```shell
spark-submit 2-sample_qc/3-population_pca_prediction.py --assign_pops
```

Plot PCA clustering for merged dataset (1000 genomes + the dataset), and for the dataset only.
You can specify the number of PCA components you want in the config file.

```shell
spark-submit 2-sample_qc/3-population_pca_prediction.py --pca-plot-assigned
```

4. Outliers identification

Now that we have the predicted populations that each sample belongs to we can run sample QC stratified by population and identify outliers within each population for each metric tested:
- number of SNPs
- number of deletions and insertions, insertion/deletion ratе
- heterozygosity rate, heterozygous/homozygous ratio
- number of transitions and transversions, transition/transversion ratio.

```shell
spark-submit 2-sample_qc/4-find_population_outliers.py
```

5. Filter data to exclude samples which fail QC.

The final step in sample QC is filtering the data to remove samples which are identified as failing in the previous script. <!At this stage samples failing on FREEMIX score and on identity checks are also removed. This samples should be in files in the annotations directory: `verify_bam_id_result_concat.selfSM` lists sample ID and FREEMIX score and `sanger_samples_excluded_after_gtcheck.txt` lists samples failing identity checks. If no samples fail identify checks the latter file could be empty.> These samples are saved in `samples_failing_qc.tsv.bgz` in the annotation directory.

```shell
spark-submit 2-sample_qc/5-filter_fail_sample_qc.py
```

### 3. Variant QC

Variant QC requires a ped file detailing any trios. This is an unheaded, tab-delimited file that contains the following columns:
- Family ID
- Proband ID
- Paternal ID
- Maternal ID
- Proband sex (1 or 2)
- Proband affected status (0 or 1)


The first step of variant QC produces a truth set of known true positive variants and generates annotation.

```shell
spark-submit 3-variant_qc/1-generate_truth_sets.py --all
```

Next an input table is generated to run the random forest on.

```shell
spark-submit 3-variant_qc/2-create_rf_ht.py
```

Next train the random forest model, note the use of --master local[*] here which is needed to ensure that the worker nodes pick up the correct python modules.

```shell
spark-submit  3-variant_qc/3-train_rf.py
```

The random forest model ID (called runhash previously, so you can find this term in the code)
will be printed to STDOUT.
It is an 8-character string consisting of letters and numbers.
Put this ID in the config file in the `rf_model_id:` parameter under the `general` section.

**Note:**
In old gnomAD releases, some random forest related functions,
(for example `train_rf_model()`) were a bit buggy
and could work incorrectly in the parallel SPARK environment.
If VariantQC functions fail with some weird message
(no space left on the device, wrong imports, etc),
try running model training on the master node only by adding `--master local[*]`
to the `spark-submit` parameters.

Now apply the random forest to the entire dataset.

```shell
spark-submit 3-variant_qc/4-apply_rf.py
```

Annotate the random forest output with metrics including synonymous variants, family annotation, transmitted/untransmitted singletons and gnomAD allele frequency. Synonymous variants are required in a file generated from VEP annotation and in the following format:

```
chr10   100202145   rs200461553 T   G   synonymous_variant
chr10   100204510   rs2862988   C   T   synonymous_variant
chr10   100204528   rs374991603 G   A   synonymous_variant
chr10   100204555   rs17880383  G   A   synonymous_variant
```

```shell
spark-submit 3-variant_qc/5-annotate_ht_after_rf.py
```

Add ranks to variants based on random forest score, and bin the variants based on this.

```shell
spark-submit 3-variant_qc/6-rank_and_bin.py
```

Create plots of the binned random forest output to use in the selection of thresholds. Separate thresholds are used for SNPs and indels.

```shell
spark-submit 3-variant_qc/7-plot_rf_output.py
```

After examining the plots and selecting suitable thresholds for SNPs and indels you can calculate the number of true positive and false positive variants remaining at your chosen thresholds and at the bins surrounding it using the following (where snp_bin and indel_bin are the thresholds selected for SNPs and indels respectively).

```shell
spark-submit 3-variant_qc/8-select_thresholds.py --snv snp_bin --indel indel_bin
```

Filter the variants in the Hail MatrixTable based on the selected threshold for SNPs and indels

```shell
spark-submit 3-variant_qc/9-filter_mt_after_variant_qc.py --snv snp_bin --indel indel_bin
```

### 4. Genotype QC

Genotype QC may be performed using a range of filters defined in *config/inputs.yaml* or it may be performed using a single set of filters for DP, GQ, VAF and bins. Here we apply a range of filters (relaxed, medium and stringent) which we have defined in the *config/inputs.yaml* file

```shell
spark-submit 4-genotype_qc/1a-apply_range_of_hard_filters.py
```

In order to evaluate these filters, variant counts per consequence per sample and transmitted/untransmitted ratio of synonymous singletons (if trios are present in the data) are calculated as follows. VEP annotation is required for this step in the following format:

```
chr10   100199947   rs367984062 A   C   intron_variant
chr10   100199976   rs774723210 G   A   missense_variant
chr10   100200004   .   C   A   missense_variant
chr10   100200012   rs144642900 C   T   missense_variant
chr1    100200019   .  C    A   stop_gained&splice_region_variant
```

```shell
spark-submit 4-genotype_qc/2-counts_per_sample.py
```

Export the filtered variants to VCF.

```shell
spark-submit 4-genotype_qc/3a-export_vcfs_range_of_hard_filters.py
```

Alternatively, to export VCFs wiht only stringent hard filter, use the 3b version of the script:

```shell
spark-submit 4-genotype_qc/3b-export_vcfs_stingent_filters.py
```
