# Getting Started With WES QC Using Hail

This guide covers WES QC using Hail. It is important to note that every dataset is different and that for best results it is not advisable to view this guide as a recipe for QC.
Each dataset will require careful tailoring and evaluation of the QC for best results.

## Before you start

In order to run through this guide you will need an OpenStack cluster with Hail and Spark installed.
It is recommended that you use `osdataproc` to create it.
Follow the [Hail on SPARK](hail-on-spark.md) guide to create such a cluster.

The ability to run WEQ-QC code on a local machine is under development.

This guide also requires a WES dataset joint called with GATK and saved as a set of multi-sample VCFs.
If starting with a Hail matrixtable, then start at [Step 2](#2-sample-qc).

## Set up

Clone the repository using:
```shell
git clone https://github.com/wtsi-hgi/wes-qc.git
cd wes_qc
```

If you are running the code on a local machine (not on the Hail cluster),
set up virtual environment using `uv`.

```bash
pip install uv # Install uv using your default Python interpreter
uv sync # install all required packages
```

Activate your virtual environment
```bash
source .venv/bin/activate
```

**Note**: Alternatively, you can work without activated virtual environment.
In this case you need to use `uv run` for each command.
For example, to run tests: `uv run make integration-test`.

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
All other paths are specified as relatives, so you won't need to edit it.

If you have to make any dataset-specific operations, you can create your own branch
and add all the code you need to it.

## How to run the code

Start a new tmux session and modify environment variables
to include the directory you originally cloned the git repo into.
To get a correct Python path, you need to have the virtual environment activated.

```shell
export PYTHONPATH=$PYTHONPATH:$(pwd)
export PYSPARK_PYTHON=$(which python)
export PYSPARK_DRIVER_PYTHON=$(which python)
```

### Manually running the code on a local machine

To manually run the code on a local machine,
run the Python and provide the path to the pipeline script:

```shell
python 1-import_data/1-import_gatk_vcfs_to_hail.py
```

### Manually run the code on a Hail cluster

To submit the jobs on a Hai cluster, you need to run the pipeline script via `spark-submit`.
For example:

```shell
spark-submit 1-import_data/1-import_gatk_vcfs_to_hail.py
```

### Automatically sync and run the code on the cluster from a local machine

If you want to modify the code on your local machine,
and then run it on the cluster, you can use two scripts
in the `scripts`

* `hlrun_local` - runs the Python script via `spark-submit`. You need to run it on the spark master node on your cluster.
* `hlrun_remote` - runs the code on the Spark cluster form your local machine.
  It performs a series of operations:
  * Sync the codebase to the remote cluster, defined by the environment variable `$hail_cluster`.
    The variable can contain the full host definition (`user@hostname`) or only hostname from the SSH config file.
  * Create tmux session on the remoter cluster
  * Run the Python script via `hlrun_local`
  * Attach to the tmux session to monitor the progress


**Warning**

The `hlrun_remote` is designed to work with only one tmux session.
To start a new task via `hlrun_remote`, first end the existing tmux session, if it exists.


## Analyze your data

### 0. Resource Preparation
All steps in this section need to be run only once before your first run. It prepares the reference dataset for the subsequent steps.

1. **Create the 1000G population prediction resource set**.

This resource set is required for the super-population prediction on the population PCA step.
Then you can reuse it with any data cohort.

```shell
spark-submit 0-resource_preparation/1-import_1kg.py --all
```

### 1. Load data

1. **Load VCFs into Hail and save as a Hail MatrixTable**

```shell
spark-submit 1-import_data/1-import_gatk_vcfs_to_hail.py
```

2. **Annotate metadata**

This script annotates samples with all provided metadata:
VerifyBamId Freemix score, self-reported sex, self-reported ethnicity, etc.

If you don't have some (or even any) of these annotations,
put `null` instead of the corresponding filename in the config file:
* `verifybamid_selfsm` - for the VerifyBamID freemix data


For each available annotation, the script prints out the list of samples that don't have annotations.
For the Freemix score it performs validation and saves the Freemix plot.

```shell
spark-submit 1-import_data/2-import_annotations.py
```

**Note:** Now, only Freemix annotation is supported. Work on other annotations is in progress.

**Note:** The output of this step is not used in further steps yet.
You can use it to run Freemix validation and get plots and sample statistics.
Using annotated data in further steps will be implemented after implementing GtCheck validation.


### 2. Sample QC

1. **Run sex imputation**

```shell
spark-submit 2-sample_qc/1-hard_filters_sex_annotation.py
```

2. **Identify samples from related individuals with PCRelate**
This step outputs a relatedness graph, a table of total statistics of relatedness and a list of related samples.
Please see config files "prune_pc_relate" for more details.

```shell
spark-submit 2-sample_qc/2-prune_related_samples.py
```

While this step identifies related samples, we keep them in the dataset since step 2.3 uses PCA score projection for population clustering. The relatedness information can be used to validate pedigree data and detect sample mislabeling.


3. **Predict populations**

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

4. Identify outliers

Now that we have the predicted populations that each sample belongs to,
we run sample QC stratified by population and identify outliers.

We test the following metrics, calculated by Hail:
* number of SNPs
* heterozygosity rate, heterozygous/homozygous ratio
* number of transitions and transversions, transition/transversion ratio.
* number of deletions and insertions, insertion/deletion ratе

For metric description, see the
[Hail sample_qc()](https://hail.is/docs/0.2/methods/genetics.html#hail.methods.sample_qc)
function description.

```shell
spark-submit 2-sample_qc/4-find_population_outliers.py
```

WES-QC pipeline identifies outliers using the gnomAD function
[`compute_stratified_metrics_filter()`](https://broadinstitute.github.io/gnomad_methods/api_reference/sample_qc/filtering.html#gnomad.sample_qc.filtering.compute_stratified_metrics_filter).
By default, this function designates as outliers any samples
that deviate more than 4 Median Absolute Deviations (MAD)
from the average by any metric.

If you need to adjust this behavior,
modify the `compute_stratified_metrics_filter_args` section in the configuration file.
Any parameters added to this section are transferred to the `compute_stratified_metrics_filter()` function.
For example, you can use the `metric_threshold` dictionary to specify individual thresholds for some metrics.

The script outputs the full list of samples with calculated metrics
(the `stratified_sample_qc`:`output_text_file` config parameter),
statistics and outlier intervals for all metrics in JSON format
(the `stratified_sample_qc`:`output_globals_json_file` config parameter),
and metric distribution plots (the `plot_sample_qc_metrics`:`plot_outfile` config parameter).

5. **Filter out samples which fail QC**

The final step in sample QC is filtering the data to remove samples which are identified as failing in the previous script.
<!At this stage samples failing on FREEMIX score and on identity checks are also removed.
This samples should be in files in the annotations directory: `verify_bam_id_result_concat.selfSM`
lists sample ID and FREEMIX score and `sanger_samples_excluded_after_gtcheck.txt` lists samples failing identity checks.
If no samples fail identify checks the latter file could be empty.
> These samples are saved in `samples_failing_qc.tsv.bgz` in the annotation directory.

```shell
spark-submit 2-sample_qc/5-filter_fail_sample_qc.py
```

### 3. Variant QC

Variant QC requires a ped file detailing any trios.
This is an unheaded, tab-delimited file that contains the following columns:
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

Next, train the random forest model.

```shell
spark-submit  3-variant_qc/3-train_rf.py
```

The random forest model ID (called _runhash_ previously, so you can find this term in the code)
will be printed to STDOUT.
It is an 8-character string consisting of letters and numbers.
Put this ID in the config file in the `rf_model_id:` parameter under the `general` section.

**Note:**
In old _gnomAD_ releases, the function `train_rf_model()`
could work incorrectly in the parallel SPARK environment.
If and VariantQC step fails with some weird message
(no space left on the device, wrong imports, etc),
try running model training on the master node only by adding `--master local[*]`
to the `spark-submit` parameters.

Now apply the random forest to the entire dataset.

```shell
spark-submit 3-variant_qc/4-apply_rf.py
```

Annotate the random forest output with metrics including synonymous variants, family annotation,
transmitted/untransmitted singletons, and gnomAD allele frequency.
Synonymous variants are required in a file generated from VEP annotation and in the following format:

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

After examining the plots and selecting suitable thresholds for SNPs and indels,
you can calculate the number of true positive and false positive variants
remaining at your chosen thresholds and at the bins surrounding it using the following
(where snv_bin and indel_bin are the thresholds selected for SNVs and indels respectively).

```shell
spark-submit 3-variant_qc/8-select_thresholds.py --snv snv_bin --indel indel_bin
```

Filter the variants in the Hail MatrixTable based on the selected threshold for SNVs and indels

```shell
spark-submit 3-variant_qc/9-filter_mt_after_variant_qc.py --snv snv_bin --indel indel_bin
```

### 4. Genotype QC

Genotype QC are performed using a range of filters defined in `config/inputs.yaml`.

To perform genotype QC, you need to determine the best combination of hard filters,
to save "good" variations as much as possible,
and get rid of all "bad" variants and genotypes at the same time.

The first script of the genotype QC helps you to analyze different combinations of hard filters
and choose optimal values.

Based on the results of the VariantQC step populate the provisional values
for the SNV and indel random forest bins in the `evaluation` part of the config file.
For example:

```yaml
    snp_bins: [ 60, 75, 90 ]
    indel_bins: [ 25,50,75 ]
    gq_vals: [ 10, 15 ]
    dp_vals: [ 5, 10 ]
    ab_vals: [ 0.2, 0.3 ]
    missing_vals: [0.0, 0.5]
```

The values for the Genotype quality (gq), read depth (dp), allele balance (ab),
and missingness (minimal percentage of genotypes where this variation is defined)
are typical, so you can start with the provided values.

If your dataset contains a control sample with known high-confident variations
(usually one of the [GIAB](https://www.nist.gov/programs-projects/genome-bottle) samples),
you can use it to calculate precision/recall values.
Add the sample control name, the corresponding VCF file, and the VEP annotation to the config:

```yaml
    giab_vcf: '{resdir}/HG001_GRCh38_benchmark.interval.illumina.vcf.gz'
    giab_cqfile: '{resdir}/all.interval.illumina.vep.info.txt'
    giab_sample: 'NA12878.alt_bwamem_GRCh38DH.20120826.exome'
```

Run the hard-filter evaluation step:

```shell
spark-submit 4-genotype_qc/1-compare_hard_filter_combinations.py --all
```

The script calculates all possible combinations of hard filters, and collects
results in the subfolder in `annotation` folder, named by the RF model ID.

Review, plot, and analyze the results.
Based on the desired balance of true positives vs. false positives
and the desired precision/recall balance, choose the three combinations of filters:
relaxed, medium, and stringent.

If needed, add more values to evaluate in the config and rerun the hard filter evaluation.

Fill in the values in the `apply_hard_filters` part of the config.

Run the Genotype QC with the chosen set of filters:

```shell
spark-submit 4-genotype_qc/2-apply_range_of_hard_filters.py
```

Export the filtered variants to VCF.
Script 3a tags all variations with the corresponding filter (relaxed, medium, stringent)
removes all variants not passing the relaxed filter, and saves the resulting data to VCF files.

```shell
spark-submit 4-genotype_qc/3a-export_vcfs_range_of_hard_filters.py
```

Alternatively, to export VCFs with only passing stringent hard filter, use the 3b version of the script:

```shell
spark-submit 4-genotype_qc/3b-export_vcfs_stingent_filters.py
```

If you want to additionally evaluate the filter statistics
(variant counts per consequence per sample and
transmitted/untransmitted ratio of synonymous singletons (if trios are present in the data))
use the following script.
VEP annotation is required for this step in the following format:

```
chr10   100199947   rs367984062 A   C   intron_variant
chr10   100199976   rs774723210 G   A   missense_variant
chr10   100200004   .   C   A   missense_variant
chr10   100200012   rs144642900 C   T   missense_variant
chr1    100200019   .  C    A   stop_gained&splice_region_variant
```

```shell
spark-submit 4-genotype_qc/4-counts_per_sample.py
```
