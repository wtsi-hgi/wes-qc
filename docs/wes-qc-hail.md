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

## How to run the code (tmux)

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


## How to run the code (jupyter notebook)
Alternatively, you can run the code in a jupyter notebook where all the steps are arranged in a sequence and divided into sections (e.g. 0-resource_preparation, 1-import_data, 2-sample_qc, 3-variant_qc, 4-genotype_qc).

The notebook is located as `scripts/run-wes-qc-pipeline-all-steps.ipynb`.

It uses hlrun_local to run the code, which will output the log file to the current directory, with the prefix of the step name, e.g. `hlrun_3-1-generate-truth-sets_20250102_125729.log`.

To run with jupyter notebook, you could either run the notebook on the cluster or on your local machine. By default, it will run on the local machine,

### Run with jupyter notebook on the local machine
Note if you are doing this on a local machine, you may need to keep your laptop on for a long time. Running on a virtual machine is recommended.

1. First, prepare a python environment with jupyter notebook installed. There is no additional dependency to install.
A quick way to do this is to install `uv` - `curl -LsSf https://astral.sh/uv/install.sh | sh`, see [uv](https://astral.sh/uv/) for more details.
Then run `uv add jupyterlab`, the `jupyter` executable will be installed in `.venv/bin/jupyter`.

2. Then, set up a SSH connection to the cluster, for example, if the cluster IP address is `172.27.1.1`, add the following to your `~/.ssh/config` file:
```
Host wes
    HostName 172.27.1.1
    User ubuntu
    IdentityFile ~/.ssh/id_rsa
```
id_rsa is the private key for accessing the cluster when the cluster was created.

3. Then, in the python environment, run the following command to start the notebook:
```
jupyter notebook scripts/run-wes-qc-pipeline-all-steps.ipynb
```
or if you are using VSCode, you can open the notebook directly in VSCode and set the python interpreter to the one with the Jupyter notebook installed.

4. Follow the instructions in the notebook to run the pipeline.
Note that the variable `jupyter_notebook_on_cluster` needs to be set to `False` in the notebook.

### Run with jupyter notebook on the cluster
1. First, prepare a python environment with jupyter notebook installed. There is no additional dependency to install.
See the above [Run with jupyter notebook on the local machine](#run-with-jupyter-notebook-on-the-local-machine) section for more details.

2. Then, start a jupyter notebook server on the cluster:
```
jupyter notebook --no-browser --port=8889 scripts/run-wes-qc-pipeline-all-steps.ipynb
```
Note you cannot use the default port 8888 because it is already used by the built-in notebook server on the master node.
The built-in notebook server on the master node is not suitable because its Spark will claim all the resources and the pipeline will not be able to run.

If you are using VSCode, you can open the notebook directly in VSCode and set the python interpreter to the system python.

3. Follow the instructions in the notebook to run the pipeline.
Note that the variable `jupyter_notebook_on_cluster` needs to be set to `True` in the notebook.

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

Specify the corresponding input file in the config for each available annotation
(follow the links to download the sample files):

* [verifybamid_selfsm:](https://wes-qc-data.cog.sanger.ac.uk/metadata/control_set_small.verify_bam_id_result.selfSM.tsv) -
  the VerifyBamID Freemix data
* [sex_metadata_file:](https://wes-qc-data.cog.sanger.ac.uk/metadata/mlwh_sample_and_sex.txt) - self-reported sex.
  A tab-separated TSV file, having at least two columns: `sample_id` and `self_reported_sex`.

If you don't have some (or even any) of these annotations,
put `null` instead of the filename in the config file.

Run the annotation script:

```shell
spark-submit 1-import_data/2-import_annotations.py
```

For each available annotation, the script prints out the list of samples that don't have annotations.
For the Freemix score it performs validation and saves the Freemix plot.

3. **Annotate and validate GtCheck results**

The good practice for clinical samples is to make independent microarray-based genotyping
together with the exome/genome sequencing.
If you have array data, you can use `bcftools gtcheck` utility to check consistency between
sequencing and microarray genotypes.

**Skipping genotype checking**: If you don't have array data, set the `wes_microarray_mapping: none`
under the `validate_gtcheck` section of the config file.
To generate the correct output matrixtable,
you need to run this script in any case, even if you don't have any array data.

To run genotype validation, you need to provide the following files in the config
(open links to obtain the sample files):

* [wes_microarray_mapping:](https://wes-qc-data.cog.sanger.ac.uk/metadata/control_set_small.microarray_mapping.tsv):
  -- the two-column tab-separated file,
  containing the expected mapping between WES and microarray samples.
  (usually, microarray studies have separated sample-preparation protocol and separate IDS)
* [microarray_ids:](https://wes-qc-data.cog.sanger.ac.uk/metadata/control_set_small.microarray_samples.txt)
  -- the list (one ID per line) of IDs, actually found in your microarray data.
  This file is expected to have the same set of IDs as in the mapping file.
  However, sometimes array ganotyping for a particular sample fails,
  and in this case it is not present in the results.
* [gtcheck_report:](https://wes-qc-data.cog.sanger.ac.uk/metadata/control_set_small.combined.gtcheck.txt)
  -- the output of the `bcftools gtcheck` command.
  You need to remove the file header and keep only data lines.

**Note:** To run `bcftools gtcheck` and generate the report,
you most probably need to convert microarray data from FAM to VCF,
and liftover it to the GRCh38 reference.
This work should be done outside WES-QC pipeline,
and is not covered by this manual.

Run the validation script:

```shell
spark-submit 1-import_data/3-validate-gtcheck.py
```

**Gtcheck validation results and interpretation**:

The validation script implements complicated logic to ensure correctness of all data.

At first, it validates the consistency of the mapping file and samples present in the data.
The script reports the IDs present in the mapping but not present in the real data, and opposite.
Also, it returns all duplicated IDs in the mapping.
After validation, the script removes from the mapping table all microarray IDs not found in the data.

Next, the script loads gtcheck table and runs a decision tree to split samples into passed and failed.

![genotype checking decision tree](pics/gtcheck_decision_tree.png)

On each decision tree step, samples are marked by the specific tag in the `validation_tags` column.

The script exports the final table for all samples, and a separate file for the samples failed validation
under the gtcheck validation dir (specified in the config file in `gtcheck_results_folder` entry).
The tags in the `validation_tags` column allows tracking the chain of decisions for each sample.
The same mechanism allows developers to extend this script and add more decision steps if needed.

Here are all already implemented tags:

| tags                                                         | Description                                                                                                              |
|--------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------|
| best_match_exist_in_mapfile, best_match_not_exist_in_mapfile | The matching gtcheck sample with the best score exists/not exists in mapping                                             |
| best_match_matched_mapfile, best_match_not_matched_mapfile   | The matching gtcheck sample with the best score is consistent/not consistent with the mapping file                       |
| score_passed,  score_failed                                  | Gtcheck score for the best matched sample passed/failed threshold check                                                  |
| mapfile_unique, mapfile_non_unique,                          | The matching array sample is unique/not unique in the mapping file                                                       |
| mapfile_pairs_have_gtcheck, no_mapfile_pairs_have_gtcheck    | There is at least one/there are no samples form the mapping file that were reported in the Gtcheck best matching samples |



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
    indel_bins: [ 25, 50, 75 ]
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

If you don't have a GIAB sample, put null in the `giab_sample` section.
The precision/recall calculations will be skipped in this case.

_Note_: For now, you still need to specify the GIAB VCF and cqfile, even if you're skipping
the precision calculation.

Run the hard-filter evaluation step:

```shell
spark-submit 4-genotype_qc/1-compare_hard_filter_combinations.py --all
```

The script calculates all possible combinations of hard filters.
Depending on the dataset size and number of evaluated combinations, the calculation can take significant time.
The script prints elapsed time and estimated time to complete after each step.

After finishing the calculations, the script
makes interactive plots in the `plot` directory with the `hard_filter_evaluation` prefix.
Also, the script saves results in the subfolder in `annotation` folder, named by the RF model ID.

**At this step, you MUST review and analyze the results to choose correct values for hardfilter combinations**.
The values for the public datasets are not suitable for your data.

All plots are interactive, and you can use the following options to explode your data:
* Zoom in/out
* Move to a specific location by dragging a graph content.
* Use Checkboxes to filter data by DP, GQ AB, and call_rate hardfilters.
  This option is especially useful when you select/deselect a checkbox and observe how your data are changed.
* Use sliders to filter by minimum/maximum bin value.
* Change between several available color maps using the dropdown menu.

If you need to analyze more data points, add required values to the config file and rerun the evaluation.
The evaluation script dumps intermediate results for each filter combination and calculates only new
combinations, so the real calculation time will be smaller than estimated.

If you need to recalculate some combinations, go to the folder with dumped results (`json_dump_folder`)
and manually delete all combinations that you need to re-evaluate.
To rerun all calculations from scratch, delete the dump folder entirely.


Based on the desired balance of true positives vs. false positives
and the desired precision/recall balance, choose the three combinations of filters:
relaxed, medium, and stringent.
Fill in the values in the `apply_hard_filters` part of the config.
If needed, add more values to evaluate in the config and rerun the hard filter evaluation.

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
