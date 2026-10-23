# Exome and Genome cohorts QC using WxS-QC pipeline

This guide explains how to perform QC of human exome/genome data cohorts
using the WxS-QC pipeline.

The general information regarding the pipeline capabilities and limitations
is available in the [readme](../README.md).

## Pre-requirements

Input data should be a multi-sample VCF file or a set of VCF files.
We tested the pipeline on variants calling by GATK4 suite using haplotype caller and joint calling steps.

The pipeline has no direct dependencies of GATK,
but it uses variant-level and genotype-level metrics calculated by GATK.
By altering the inputs to the RF model, the variant QC could be adapted for other variant callers
(FreeBayes, Strelka2, etc).
Details are provided in the [Variant QC section](#stage-3-variant-qc).

The adaptation of the variant QC for modern neural network-based callers,
like DeepVariant and DRAGEN, is underway.

We suggest not applying the GATK4 VQSR,
because the pipeline has its own more flexible variant QC block.

### Public example dataset

We use the open set of data from 1000 genomes project as an example of data analysis.
You can download and use these data to ensure that your installation is running properly.

**It is important to note that every dataset is different and
requires careful tailoring and evaluation of the results.**
For the best results, don't consider this guide as an end-to-end recipe for QC.
Review the results and tune the analysis parameters when we suggest that you do it.

## Setting up

### Set up pipeline code and environment

To run through this guide, you need a machine or a cluster with Python and Hail installed.
To set up the environment, refer to the [setup howto](wxs-qc_setup.md).

**Note:**
If you have Hail installed on a SPARK cluster,
use `spartk-submit` instead of regular `pyhton`
to run pipeline scripts.

### Obtain resource files

The WxS-QC pipeline uses a set of resource data.
For the resource description and obtaining instruction,
refer to the [resources description section](wxs-qc_prepare-resources.md).

### Obtain public example dataset

We prepared the open set of data from the 1000 genomes project as an example of data analysis.
You can download and use these data to ensure that your pipeline installation is running properly.
You can download the dataset from the WSI page:

```bash
wget https://wxs-qc-data.cog.sanger.ac.uk/wxs-qc_public_dataset_v3.tar
tar xvf wxs-qc_public_dataset_v3.tar
```

**It is important to note that every dataset is different and
requires careful tailoring and evaluation of the results.**
For the best results, don't consider this guide as an end-to-end recipe for QC.
Review the results and tune the analysis parameters when we suggest that you do it.


## Prepare the data and configuration

The WxS-QC config file is a YAML file with the ability to reference one field from another.
Each string in brackets like `{mtdir}` or `{general.dataset_name}` is a reference to another config section.
See the config file caption for details.

### Make the configuration file for your dataset

Create a new config file for your dataset.

By default, all scripts will use the config file named `inputs.yaml`.
You can make a symlink for it to keep the config name meaningful.

```shell
cd config
cp public-dataset.yaml my_project.yaml
ln -snf my_project.yaml inputs.yaml
cd ..
```

If you want to test your WxS-QC installation with the open test data,
you need to modify only the `data_root` parameter.

To start processing your own data, check and modify the following parameters:

* In the `general` section
  * `data_root` and `dataset_name` specifying the name of your dataset and path to it.
    Give a meaningful name to your dataset.
    By default, `data_root` variable contains reference to the dataset name, but you can use an exact file path if you want.
  * `tmp_dir` by default references to the local folder `tmp` inside your dataset folder.
    You can change
  * `onekg_resource_dir`: The place for the 1000-Genome VCFs.
    By default, it is `resources/mini_1000G` field under your data root.
    If you work with your own data, change it to downloaded full-sized 1000 genomes data, as described in
    the [resources howto](wxs-qc_prepare-resources.md).
  * `rf_model_id` leave it empty for now and specify after creating the random forest model
    on the Variant QC stage
* `stage0 -> indir` and `stage0 -> kg_pop_file` the directory for the 1000G genomes sample data. See below how to obtain it.

All other files and resources you need are described in the corresponding sections
of this manual.

### Create the analysis folder for your data

The default settings in the config file assume that all your data and analysis results
will be stored in a specific analysis folder.

To create the folder for your dataset with all required subfolders, you can run the script:

```shell
python stage0_resource_preparation/res0_create_project_folder.py
```

The script will take all values from the config file and create the dataset folder
with the following subfolders:
* `annotations` - various table reports generated by the pipeline
* `matrixtables` - Hail matrixtables
* `metadata` - all metadata for your dataset: self-reported sex, self-reported ethnicity, etc
* `plots` - plots generated by the pipeline
* `resources` - set of resource data
* `tmp` - temporary data folder
* `variant_qc_random_forest`
* `vcf_afterqc_export` - for exporting final VCFs after QC

Copy/symlink to the folder your resources and data.

## Stage 0. Prepare resource data

Steps in this section prepare the reference dataset for the subsequent steps.
You only need to run it once.
For next analysis you can symlink resource matrixtables.

### Create the 1000G population prediction resource set

This resource set is required for the super-population prediction on the population PCA step.
It converts the 1000 genomes VCFs into the matrixtable.

```shell
python stage0_resource_preparation/res1_import_1kg.py --all
```

### Create the combined Truth Set table

Run this step to combine all available variation resources (1000 Genomes, Mills, Hapmap, etc)
into a single table of truth variants.

```shell
python stage0_resource_preparation/res2_generate_truthset_ht.py
```

### (Optional) Extract whole-genome gnomAD frequencies

The Variant QC part of the pipeline uses variant frequencies from the
[gnomAD project](https://gnomad.broadinstitute.org/)
to find _de novo_ variations.
The resource bundle contains variant frequencies for exome variants.
If you want to use whole-genome frequencies, refer to the
[resource preparation](wxs-qc_prepare-resources.md#using-whole-genome-gnomad-variant-frequencies) section
to prepare whole-genome gnomAD frequencies.

## Stage 1. Load data

This stage focuses on the variations loading and validation against external data.

### Load VCFs into Hail and save as a Hail MatrixTable

Specify in the config file under the `stage1 -> gatk_vcf_indir` data entry
the path to the directory that you created for pre-QC VCFs.

Run data import:

```shell
python stage1_import_data/imt1_import_vcfs.py
```

### Annotate metadata

This script annotates samples with the provided metadata:
VerifyBamId Freemix score and self-reported sex.

Specify the corresponding input file in the `general->metadata` config for each available annotation
(follow the links to download the sample files):

* [verifybamid_selfsm—:](https://wes-qc-data.cog.sanger.ac.uk/metadata/control_set_small.verify_bam_id_result.selfSM.tsv) - the VerifyBamID Freemix data.
  To prepare this file, join together results of the individual VerifyBamID runs.
* [sex_metadata_file:](https://wes-qc-data.cog.sanger.ac.uk/metadata/mlwh_sample_and_sex.txt) - self-reported sex.
  A tab-separated TSV file, having at least two columns: `sample_id` and `self_reported_sex`.
  The `sample_id` column contains IDs of your samples (same as in your input VCFs).
  The `self_reported_sex` contains sex definition: `female`, `male` or `undefined`.
* _batch_ - a two-column tab delimited file specifying batch for each sample.
  If present, WxS-QC applies sample QC statistics for each batch individually.
  Is useful when your data contain subsets sequenced using different wel-lab protocols,
  The example dataset contains no batches.

If you don't have some (or even any) of these annotations,
put `null` instead of the filename in the config file.

You can find examples of metadata files
in the `tests/test_data/metadata` folder created on the [Obtain resource files](#obtain-resource-files) step.

Run the annotation script:

```shell
python stage1_import_data/imt2_import_annotations.py
```

For each available annotation, the script prints out the list of samples that don't have annotations.
For the Freemix score it performs validation and saves the Freemix plot.

The `remove_freemix_outliers` option in the `validate_verifybamid` config section
controls whether we remove all Freemix outliers at this step
or pass it to the sample QC.

### Annotate and validate GtCheck results

A good practice for clinical samples is to make independent genotyping
together with the exome/genome sequencing, for example, using SNP microarrays.
If you have any independent genotype data, you can run the `bcftools gtcheck` utility
and validate consistency between genotypes using the WxS-QC pipeline.
To do this, refer to the [Validate GtCheck results](wxs-qc_gtcheck-validation.md) howto.

The external genotype validation step is optional.
However, to generate the correct output,
you need to run this script in any case, even if you don't have any external genotype data.
To do it, set the `wes_microarray_mapping: none`
under the `validate_gtcheck` section of the config file, and run the step:

```shell
python stage1_import_data/imt3_validate_gtcheck.py
```


### Plot mutations spectra

Plotting the mutation spectra can help you to identify batch-level artifacts.
To do it, run the calculation script:

```shell
python stage1_import_data/imt_4_mutation_spectra_preqc.py
```

The script saves the plot in the HTML file specified under the
`plot_mutation_spectra_preqc`:`mut_spectra_path` config section.
Also, you can specify the IQR range for outliers and change the plot size if needed.


## Stage 2. Sample QC

The WxS-QC pipeline implements three different methods for sample QC.
You can choose it in the `general -> sample_qc_method` section of the config file:
  * `pop` - Stratified QC with PCA-assigned superpopulations with PCA projection (from gnomAD v3)
  * `nn` - Nearest neighbors (from gnomAD v4)
  * `lr` - Linear regression (from gnomAD v4)

For all methods we calculate a set of key quality metrics:
* number of SNPs
* heterozygosity rate, heterozygous/homozygous ratio
* number of transitions and transversions, transition/transversion ratio.
* number of deletions and insertions, insertion/deletion rate

We flag as outliers the samples that deviate significantly (typically more than four Median Absolute Deviations)
from the median of their reference group.
The MAD level can be configured individually per metric on the step
[Identify outliers](#identify-outliers).
The difference between methods is how they define a reference group for each sample.
Depending on the chosen method, the step outputs can slightly differ.
You can try all of it and compare the results.

The control samples, named in the `general -> metadata -> control_samples` section,
are kept in the dataset despite any results of the sample QC check.

**Warning:**
All sampleQC methods are designed to work with relatively big population-level datasets.
We recommend having >=120 independent samples (200+ per ancestry/stratum is preferable) to obtain reliable results.

#### Superpopulation Stratified Outlier Detection (pop)
This method, originally used in gnomAD v3, performs outlier detection within predefined superpopulation groups
(e.g., EUR, AFR, EAS, ...), identified on the
PCA-based super-population prediction.

This approach works the best for groups with a similar genetic background,
which fits into a reference superpopulations from 1000 Genomes data.
More information can be found in the [gnomAD v3 documentation](https://gnomad.broadinstitute.org/news/2019-10-gnomad-v3-0/).
For metric description, see the
[Hail sample_qc()](https://hail.is/docs/0.2/methods/genetics.html#hail.methods.sample_qc)
function description.
You can set up the number of PCA components for the script
`3-population_pca_prediction.py`
in the section `stage2 -> pop_pca` of the config file.

**Warning:** The PCA superpopulation prediction and stratified filtering can incorrectly mark as outliers
samples whose ancestries are not represented in the 1000 Genomes data
(which is a common case for non-European populations).
Also, stratified filtering produces false outliers for superpopulations with <100 samples in each subpopulation.
Use the _pop_ method only for large cohorts with ancestries covered by 1000 Genomes reference panel.
The _Nearest Neighbors_ and _Linear Regression_ methods were specially introduced
to correctly perform sample QC for samples with ancestries not represented in the panel.

#### Nearest Neighbors (nn)
The Nearest Neighbors method, introduced in gnomAD v4,
offers a more granular approach to population stratification.
Instead of grouping samples into pre-defined superpopulations,
it identifies a set of "genomic neighbors" for each sample based on their coordinates in Principal Component (PC) space.
Then it compares QC metrics for a given sample against the distribution of those same metrics among its closest genetic neighbors.

This method is particularly effective for samples with admixed or rare ancestries
that do not fit neatly into major superpopulations.
Further details on this approach are available in the [gnomAD v4 post](https://gnomad.broadinstitute.org/news/2023-11-gnomad-v4-0/).

You can set ll the parameters for this function in the section
`stage2 -> stratified_sample_qc -> determine_nearest_neighbors_args`


#### Linear Regression Residuals (lr)
This method uses a regression-based approach to account for the continuous relationship between genetic ancestry
(represented by PCs) and QC metrics.
It fits a linear model to predict each QC metric from the sample's principal components as independent variables.
The "residuals" (difference between the observed metric value and the value predicted by the model)
represent the portion of the metric variation that cannot be explained by ancestry.
Samples with extreme residuals are flagged as outliers.

This method provides a robust way to normalize metrics across a continuous ancestry spectrum.
Relevant technical details can be found in the
[gnomAD methods documentation](https://broadinstitute.github.io/gnomad_methods/api_reference/sample_qc/filtering.html#gnomad.sample_qc.filtering.compute_qc_metrics_residuals).

### Run sex imputation

```shell
python stage2_sample_qc/sqc1_sex_annotation.py
```

The script imputes genetic sex for all samples, and saves the results in the
`sex_annotated.sex_check.tsv.bgz` table in the annotation folder.
Then the script identifies F-stat outliers, and saves it in the
`sex_annotation_f_stat_outliers.tsv` table.
Finally, if self-reported sex is available, the script identifies samples that have
a conflict between self-reported sex and genetically imputed sex, and saves it in the table
`conflicting_sex.tsv`.
The `undefined` sex is treated as confliscting both for self-reported and for imputed sex.


### Identify samples from related individuals with PC-Relate

This script identifies related individuals in the dataset
to ensure they do not bias population structure analysis in later steps.

To identify related samples
and to perform stratified outlier detection on the subsequent steps,
we use a subset of high-quality variants:

* SNVs, autosomal, and biallelic.
* called in >= 99% of the samples (configured by `filter_params -> call_rate_threshold` parameter),
* with >= 5% of alternative alleles (configured by `filter_params -> af_threshold` parameter),
* with Hardy-Weinberg Equilibrium (HWE) p_value >= 1E-5 (controlled by `filter_params -> hwe_threshold` parameter).
  The HWE p-value is calculated using the
  [Hail VariantQC function](https://hail.is/docs/0.2/methods/genetics.html#hail.methods.variant_qc).

The challenge is that both relatedness estimation and population structure analysis
can interfere with each other: related samples distort PCA, and population structure
can cause false relatedness signals.
To address this, we adopted the approach described in
[SISG Workshop](https://uw-gac.github.io/topmed_workshop_2018/pc-relate.html).

The script uses a two-pass approach:

1. **KING kinship estimation**: First, the script runs LD pruning and uses the
   [KING](https://hail.is/docs/0.2/methods/relatedness.html#hail.methods.king)
   kinship estimator to get an initial split of samples into related and unrelated sets.
   KING works reliably on any dataset.
2. **PCA on unrelated samples + projection to related**:
   PCA is computed on the unrelated set, and related samples are projected onto this PCA space.
3. **PC-Relate**: Using the PCA scores from the previous step, the script runs
   [PC-Relate](https://hail.is/docs/0.2/methods/relatedness.html#hail.methods.pc_relate)
   to refine kinship coefficients and IBD (Identity By Descent) estimates.
   PC-Relate uses PCA to account for population structure, which improves
   relatedness detection across different superpopulations.
4. **Final split**: Based on the refined PC-Relate scores, the script identifies
   a **maximal independent set** of unrelated individuals and generates
   a final list of related samples that exceed the kinship threshold.
5. **Final PCA**: PCA is computed on the final unrelated set,
   and related samples are projected onto this PCA space.
   These PCA scores are used on the subsequent steps to perform stratified sample filtering

This two-pass approach correctly identifies related samples
even in datasets with diverse superpopulations.

> **Under validation: relatedness on small cohorts.**
> This approach was validated on large cohorts, where it behaves as described above.
> Its behaviour on small cohorts is still under validation and should not be relied on.
> On the small trio test dataset (`control_set_small_v2`, 14 samples) PC-Relate reports
> kinship around 0.10 for known parent-offspring pairs instead of the expected ~0.25.
> No pair then clears the default `relatedness_threshold` of 0.125, so the final list of
> related samples comes out empty and the final PCA treats every sample as unrelated.
> PC-Relate needs enough samples and variants to estimate the ancestry-adjusted
> kinship reliably, so a small cohort is the suspected cause, but this is not confirmed.
> If you run this step on a small cohort, inspect the kinship values in the relatedness
> table and the Kinship vs. IBD2 plot before trusting the related-samples list, and
> consider lowering `stage2 -> pc_relate_params -> relatedness_threshold` to match the
> kinship scale you actually observe.
> The integration-test assertions covering this behaviour are commented out for the same
> reason; see `tests/integration_tests/test_integration_trios.py` and
> [the development guide](wxs-qc_development.md).

The script is split into resumable substeps. Run the substeps in order:

```shell
python stage2_sample_qc/sqc2_identify_related_samples.py --filter-mt
python stage2_sample_qc/sqc2_identify_related_samples.py --pc-relate
python stage2_sample_qc/sqc2_identify_related_samples.py --plot-pca
python stage2_sample_qc/sqc2_identify_related_samples.py --validate-trios
```

The same substeps can also be run in a single invocation:

```shell
python stage2_sample_qc/sqc2_identify_related_samples.py --all
```

The `--filter-mt` substep removes control samples and filters to high-quality variants,
The `--pc-relate` substep splits samples into KING-related and KING-unrelated sets,
runs PCA on the KING-unrelated samples, projects the KING-related samples, runs PC-Relate,
and writes the final table of related samples.
The `--plot-pca` substep plots relatedness, computes the final PCA after PC-Relate pruning,
and writes the PCA matrix, scores, and loadings used by later sample-QC steps.

The step outputs the list of related samples and PCA scores.
The relatedness information can be used to validate pedigree data and detect sample mislabeling.
It also plots the Koch's relatedness plots (Kinship vs. IBD2) and PC1/PC2 sample scatterplot.

The experimental`--validate-trios` substep cross-checks the pedigree supplied in
`general -> metadata -> pedigree` against the PC-Relate relatedness table,
and writes a revised pedigree containing only the trios it can confirm.
For each parent-child pair declared in the pedigree, it looks up the corresponding
kinship and IBD statistics in the relatedness table and keeps the pair only if it
satisfies conservative parent-child thresholds
(`kin` between 0.19 and 0.4, `ibd0` < 0.025, `ibd1` > 0.75, `ibd2` < 0.125, see
`PEDIGREE_THRESHOLDS` in `wxs_qc/constants.py`). A trio is retained in the revised
pedigree only if both the maternal and paternal pairs are confirmed this way.
This substep is skipped if `stage2 -> relatedness_output -> revised_pedigree` is set
to `null` in the config. The revised pedigree file can then be used in place of the
original pedigree in later stages (e.g. `stage3 -> pedfile`).

### Predict super-populations using PCA

This script predicts the super-population of each sample
by comparing the study dataset with a reference 1000 Genomes panel.

We perform ld_prune on the 1000 Genomes reference panel,
and them merge the manel with the dataset.
To merge samples, we intersect variants between the two datasets.
Therefore, we run PCA on the subset of high-quality variants from the experimental dataset,
and the 1000 Genomes dataset.
This ensures that the LD structure is defined by the reference,
making results more stable and comparable across studies.

```shell
python stage2_sample_qc/sqc3_pca_population_prediction.py --merge-and-ldprune
```

Next, we run PAC only for 1000 genomes samples,
and then use PC projection to predict super-populations in the experimental dataset.
This approach, compared to the direct PCA of the merged dataset,
makes PC axes more stable and not distorted by related samples or unusual ancestry.
As a result, we can handle datasets with any number of related individuals,
and PCs from different studies are more comparable.

```shell
python stage2_sample_qc/sqc3_pca_population_prediction.py --pca
```

Plot 1KG PCA. In this substep, all dataset samples should be labelled as `N/A`.

Since the set of variations for the PCA run is unique for each dataset,
the PCA results for 1000 genomes can differ between runs. However, for the suggested cohort sizes,
the number of high-quality variants in the dataset is big enough to make PCA results visually comparable.

```shell
python stage2_sample_qc/sqc3_pca_population_prediction.py --pca-plot
```
In this substep, all dataset samples should be labelled as `N/A`.

Project the study samples onto these principal components,
run a **Random Forest classifier** to assign populations,
and plot the results.
The plotting substep shows PCA clustering for the merged dataset (1000 Genomes + the dataset) and for the dataset only.
You can specify the number of PCA components you want in the config file.

The projection ensures that the study cohort’s structure is measured against
a known global reference without biasing the principal components themselves.

```shell
python stage2_sample_qc/sqc3_pca_population_prediction.py --assign_pops
python stage2_sample_qc/sqc3_pca_population_prediction.py --pca-plot-assigned
```
The script outputs population assignment tables and PCA plots.
Review graphs to ensure that the superpopulation structure matches your expectations.

### Identify outliers

Now we run sample QC stratified by the reference group (depending on the chosen method) and identify outliers.

We test the following metrics, calculated by Hail:
* number of SNPs
* heterozygosity rate, heterozygous/homozygous ratio
* number of transitions and transversions, transition/transversion ratio.
* number of deletions and insertions, insertion/deletion ratе

For metric description, see the
[Hail sample_qc()](https://hail.is/docs/0.2/methods/genetics.html#hail.methods.sample_qc)
function description.

```shell
python stage2_sample_qc/sqc4_find_outliers.py
```

WxS-QC pipeline identifies outliers using the gnomAD function
[`compute_stratified_metrics_filter()`](https://broadinstitute.github.io/gnomad_methods/api_reference/sample_qc/filtering.html#gnomad.sample_qc.filtering.compute_stratified_metrics_filter).
By default, this function designates as outliers any samples
that deviate more than four Median Absolute Deviations (MAD)
from the average by any metric.
If you need to adjust this behavior,
modify the `compute_stratified_metrics_filter_args` section in the configuration file.
Any parameters added to this section are transferred to the `compute_stratified_metrics_filter()` function.
For example, you can use the `metric_threshold` dictionary to specify individual thresholds for some metrics.

The script outputs the full list of samples with calculated metrics
(the `stratified_sample_qc`:`output_text_file` config parameter),
statistics, and outlier intervals for all metrics in JSON format
(the `stratified_sample_qc`:`output_globals_json_file` config parameter).
The script also plots distribution histograms for all metrics, and
saves them in the folder defined by the `plot_sample_qc_metrics`:`plot_outdir`
config parameter (a set of individual plots and one combined plot for all metrics and populations).
To change default number of bins, use the `n_bins` config parameter.

The sampleQC metrics and results for each sample are saved in the
`annotations/stratified_sample_qc.{method}.tsv.gz` table.

Depending on the chosen method, the step outputs can slightly differ:
* `pop` plots histograms for each metric and each superpopulation individually
* `nn` and `lr` plot only one histogram per metric, because each sample is compared to its own average
  Also, they store method-specific sample information in files
  `annotations/nearest_neighbours.tsv.gz` and `annotations/residuals.tsv.gz` respectively.
* `lr` additionally stores the linear regression parameters in the `annotations/lms.json` file.

### Filter out samples which fail QC

The final step in sample QC is filtering the data to remove samples which are identified as failing in the previous script.

At this step, you can provide an additional list of samples to remove.
It could be the samples failing FreeMix score, F-start check, or any other.
To use it, specify the file in the `samples_to_remove_file` section of the config.
This is a plain text file without a header, containing one sample per line.

If you don't want t remove any additional samples, put `null` in the config instead of sample file name.

**Note:** more convenient combined sample statistics are in the implementation.

```shell
python stage2_sample_qc/sqc5_filter_fail_samples.py
```

Samples that fail sample QC are saved in the file
`annotations/samples_failing_qc.tsv`

## Stage 3. Variant QC

The Variant QC stage trains and runs a random forest model to estimate variation quality
and rank all variations by this estimation.

Historically, this stage was designed for “classic” variant callers,
which produce a large set of variant-level statistics.
In the development version, we’re updating this part to support DeepVariant and other neural network-based callers,
which provide a very brief set of statistics.

To train the predicting model, we need a set of True-Positive (TP) and False-positive (FP) variations.
Because for a new dataset we don't have the real TP and FP, we use the following approach:

* For likely-true-positive variants, we use all variations that were found and reported
  in public databases: 1000 genomes, HapMap, Mills, etc (see the full list of sources on the step 0.2)
* For likely-false-positive variants, we use variants with low quality scores.
  At this point we use quality metrics provided py the variant caller:
  * [QualByDepth (QD)](https://gatk.broadinstitute.org/hc/en-us/articles/30331972014747-QualByDepth) — the QUAL score normalized by allele depth (AD) for a variant
  * [FisherStrand (FS)](https://gatk.broadinstitute.org/hc/en-us/articles/30332004699931-FisherStrand) — strand bias estimated using Fisher's Exact Test
  * [RMSMappingQuality (MQ)](https://gatk.broadinstitute.org/hc/en-us/articles/30332003370267-RMSMappingQuality) — root mean square of the mapping quality of reads across all samples

In this documentation, the TP and FP abbreviations stand for likely-true-positive and likely-false-positive.
Detailed information regarding all variant metrics produced by GATK is available here:
https://gatk.broadinstitute.org/hc/en-us/articles/30331989211419--Tool-Documentation-Index

Variant QC uses the optional pedigree file detailing trios existing in the dataset.
The file should follow the [FAM file](https://www.cog-genomics.org/plink/1.9/formats#fam)
notation from PLINK utility.
This is an unheaded, tab-delimited file that contains the following columns:
- Family ID
- Proband ID
- Paternal ID
- Maternal ID
- Proband sex (1-male, 2-female, 0-unknown)
- Proband affected status (0 or 1)

If you have it, specify the pedigree file in the `general->metadata->pedigree` section of the config file.
If you don't have pedigree data, several substeps will be skipped, and some metrics
for the final graphs won't be calculated.

### Generate family statistics

The first step of variant QC is to split multi-allelic variants and annotate it with family statistics.

```shell
python stage3_variant_qc/vqc1_split_and_family_annotate.py --all
```

### Generate RF test data
Next, we combine all the data and generate the input table to run the random forest on.

The random model uses a set of variant-level statistics, populated after the joint calling step.
We have tested the pipeline only on data produced by [GATK v4](https://gatk.broadinstitute.org/hc/en-us) caller.
Therefore, the list of statistics we use for the random foresm model is hardcoded in the file
[constants.py](../utils/constants.py)

To adapt the variant QC for your variant caller (FreeBayes, Strelka2, Octopus, etc),
put the annoations available in the `INFO` field of your VCF file in the
`INFO_FEATURES` of the `constants.py` file.

```shell
python stage3_variant_qc/vqc2_create_rf_ht.py
```

### Train the random forest model

Then, we train the RF model on the constructed dataset.
A random subset of the true-positive and false-positive training variants is withheld for testing the trained model;
the rest is used for training.

The size of this subset is controlled by the `stage3 -> rf_test_percentage` config parameter,
which specifies the percentage (0-100] of variants to withhold for testing,
chosen at random independently from each of the true-positive and false-positive sets.
By default (`stage3 -> train_rf -> seed: null`), a different random test subset is selected on every run;
set this to a fixed integer to select the same test subset across runs.

Here we use `--manual-model-id` parameter to
manually set the random forest model ID
(called _runhash_ previously, so you can find this term in the code)

```shell
python  stage3_variant_qc/vqc3_train_rf_model.py --manual-model-id wxs-qc_public
```

If you don't set the the RF ID, the script will
generate and print to STDOUT the random ID: 8-character string consisting of letters and numbers.
Put this ID in the config file in the `rf_model_id:` parameter under the `general` section.

**Note:**
The function `train_rf_model()`
could work incorrectly in the parallel SPARK environment,
depending on the Hail library version.
If a Variant QC step fails with some weird message
(no space left on the device, wrong imports, etc),
try running model training on the master node only by adding `--master local[*]`
to the `spark-submit` parameters,
or try to downgrade Hail to the previous version.

### Apply the random forest model

Now apply the random forest to the entire dataset.

```shell
python stage3_variant_qc/vqc4_apply_rf_model.py
```

### Annotate random forest model results

Annotate the random forest output with metrics including synonymous variants, family annotation,
transmitted/untransmitted singletons, and gnomAD allele frequency.
Synonymous variants are extracted from the VEP consequences file
(specified in `general->metadata->vep_consequences`).
The script automatically filters synonymous variants from this file.

The VEP consequences file should be in the following format
(note: the `Consequence` column is used to filter synonymous variants):

```
chr10   100202145   rs200461553 T   G   synonymous_variant  synonymous_variant  GENE1  HGNC:1234
chr10   100204510   rs2862988   C   T   missense_variant    missense_variant    GENE2  HGNC:5678
chr10   100204528   rs374991603 G   A   synonymous_variant  synonymous_variant  GENE1  HGNC:1234
```

If you don't have VEP annotation, you can set `general->metadata->vep_consequences` to null
This will disable calculation of singleton metrics.

```shell
python stage3_variant_qc/vqc5_annotate_ht_after_rf.py
```

### Group variants by ranks

At this step, we rank all variants depending on their RF score and
group it into 100 bins (bin 1 is the most quality variants.)

```shell
python stage3_variant_qc/vqc6_rank_and_bin.py
```

### Plot and analyse the results

Create plots of the binned random forest output to use in the selection of thresholds.
Separate thresholds are used for SNPs and indels.

```shell
python stage3_variant_qc/vqc7_plot_metrics.py
```

The Variant QC is always a trading between sensitivity and quality.
Examine the plots and choose the near-optimal region for RF bin,
that preserves as many TP variants as possible and at the same time eliminating most part of FP variants.
You can refer to other graphs to ensure that the chosen region also has the expected metrics.

At the genotype QC stage, we run the evaluation of different hardfilter combinations
that allow you to improve the results.
Therefore, at this point you don't need to make a final decision.
You only need to choose a provisional RF bin interval to analyze it on the next steps.

### (Optional) Checking threshold selection

If you want, you can calculate the number of true positive and false positive variants
remaining at your chosen thresholds before running genotype evaluation.

Steps 3.8 and 3.9 had been implemented before the hardfilter evaluation was introduced.
These scripts use only RF bin filtering and run faster than the full hardfilter evaluation.
You can use it to quickly estimate the percentage of variants passing bin filtering.
However, the hard filter evaluation provides more comprehensive information, so typically we don’t use these scripts.

Use for `snv_bin` and `indel_bin` the thresholds selected for SNVs and indels respectively).

```shell
python stage3_variant_qc/vqc8_select_thresholds.py --snv snv_bin --indel indel_bin
```

If you want to manually explore remaining variations,
you can filter the variants in the Hail MatrixTable based on the selected threshold for SNVs and indels.

```shell
python stage3_variant_qc/vqc9_filter_mt_manual_vqc.py --snv snv_bin --indel indel_bin
```

## Stage 4. Genotype QC

In the genotype QC stage we need to remove genotypes that are not high quality enough.
However, by removing genotypes that don't match certain filter thresholds,
we always remove some percentage of real existing genotypes.

To achieve the best result, we need to determine the best combination of hard filters,
to save "good" variations as much as possible,
and get rid of all "bad" variants and genotypes at the same time.

The first script of the genotype QC helps you to analyze different combinations of hard filters
and choose optimal values.

The first hard filter that we use, is the random forest bin,
determined in the variant QC stage.
This filter applies on the variation level, removing
all genotypes for the variation above the threshold
(for RF bin smaller values are better)

### Choose the set of hard filters for evaluation

Based on the results of the variant QC stage, populate the provisional values
for the SNV and indel random forest bins in the `evaluation` part of the config file.
For example:

```yaml
    snp_bins: [ 60, 75, 90 ]
    indel_bins: [ 25, 50, 75 ]
```

On the genotype level, we use the set of per-genotype hardfilters:
Genotype quality (gq), read depth (dp), and allele balance (ab).
Finally, we calculate **call rate** (also can be found in the code as *missingness*) —
the minimal percentage of genotypes where this variation remains defined
after applying per-genotype hard filters. This filter also applies on the variant-level.

For all these parameters, you can start with the following default values:

```yaml
    gq_vals: [ 10, 15 ]
    dp_vals: [ 5, 10 ]
    ab_vals: [ 0.2, 0.3 ]
    missing_vals: [0.0, 0.5]
```

If your dataset contains a control sample with known high-confident variations
(usually the first [GIAB](https://www.nist.gov/programs-projects/genome-bottle) sample **NA12878**),
you can use it to calculate precision/recall values.
Add the sample control name, the corresponding VCF file, and the VEP annotation to the config:

```yaml
    giab_vcf: '{resdir}/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz'
    giab_cqfile: '{resdir}/all.interval.illumina.vep.info.txt'
    giab_sample: 'NA12878'
```

Data files for the `GIAB HG001/NA12878` sample are available in the resource data
downloaded by the testing script
(see the [Obtain resource files](#obtain-resource-files) section).

_Note_: For now, you still have to specify the correct GIAB VCF and cqfile,
even if you're skipping the precision/recall calculation.

If you don't have a GIAB sample, put `null` in the `giab_sample` section.
The precision/recall calculations will be skipped in this case.

### Make the panel bed file for precision/recall evaluation
The precision-recall calculation compares the variations in the control GIAB VCF with the variations
actually present in your dataset GIAB sample.
If you run precision/recall calculations,
you should use
the intersection of high-confident regions between GIAB VCF (`HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz`)
and your sequencing panel.

Intersect it using `bedtools` to get the overlaps:

```shell
bedtools intersect -a HG001_GRCh38_1_22_v4.2.1_benchmark.bed -b your_exome_panel_intervals.bed > HG001_exome_intersected.bed
```

Place the resulting panel file to the `metadata` folder and add it to the config:

```yaml
stage4:
  evaluation:
    prec_recall_panel_bed: '${cvars.metadir}/HG001_exome_intersected.bed'
```

### Run the hard-filter evaluation step

```shell
python stage4_genotype_qc/gqc1_compare_hard_filter_combinations.py --all
```

The script calculates all possible combinations of hard filters.
Depending on the dataset size and number of evaluated combinations, the calculation can take significant time.
The script prints elapsed time and estimated time to complete after each filter combination.

### Evaluate step outputs

After finishing the calculations, the script
makes interactive plots in the `plot` directory with the `hard_filter_evaluation` prefix.
Also, the script saves results in the subfolder in `annotation` folder, named by the RF model ID.

**At this step, you MUST review and analyze the results to choose correct values for hardfilter combinations**.
The values for the public datasets are not suitable for your data.
Choosing the correct combination requires professional knowledge and could be tricky.
For more detailed explanation of this process you could review some relevant publications,
(for example https://pmc.ncbi.nlm.nih.gov/articles/PMC11747307/).

For each hardfilter combination, the script calculates the following metrics:
* Percentage of likely-true-positives (TP) and likely-false-positives (FP) variants
  (see the explanation above in the variant QC stage description).
* Precision and Recall, calculated on the GIAB sample (if it is present in the data).
  For the GIAB sample, we can calculate the real true-positive, false-positive, and false-negative variations.
  Therefore, we assume these data as more confident for choosing the optimal hardfilter evaluation
  The TP/FP values should be generally consistent with precision/recall.
  If the GIAB sample is not available, the script outputs -1.0 both for precision and recall
  * For Indels, the script additionally calculates the precision/recall values separately for in-frame and frameshift indels.
* If trios are available, the script also calculates the following values (if trios are not available, script outputs -1):
  * Transmitted/untransmitted SNV ratio — the ration of transmitted/untransmitted singleton variants across trio. Should be close to 1
  * The rate of mendelian errors: variations where the child's genotype does not match any of the possible expected combinations based on the parents' genotypes.

For these data, the script makes the following set of plots:
* FP vs TP
* recall vs precision
* transmitted/untransmitted ratio vs TP
* Mendelian error rate vs TP.

_Note_: currently, the script outputs all the graphs, even if the corresponding data are not available.
You can ignore graphs for the skipped parameters.

All plots are interactive, and you can use the following options to explode your data:
* Move to a specific location by dragging a graph content.
* Use checkboxes to filter data by `DP`, `GQ` `AB`, and `call_rate` hardfilters.
  This option is especially useful when you select/deselect a checkbox and observe how your data change.
* Use sliders to filter by minimum/maximum bin value.
* Change between several available color maps using the dropdown menu.
* Zoom in/out.
  For `FP/TP` and `precision/recall` graphs the aspect ratio is fixed to `1:1`.
  You can disable the aspect ratio fix using the checkbox under the bin selection sliders.

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

### Run the Genotype QC with the chosen set of filters

Now you can apply your custom thresholds and make variants with corresponding filters.
The pipeline can annotate all variations with VEP consequences.
To make the consequence file, extract it from VEP-annotated VCF using `bcftools`:
```bash
bcftools view dataset.vep-annotated.vcf.gz | bcftools +split-vep -s worst -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%CSQ\t%Consequence\t%SYMBOL\t%HGNC_ID\n' | bgzip > dataset.all_consequences_with_gene_and_csq.tsv.bgz
```
The example of the file can be found in the
[WSI bucket](https://wes-qc-data.cog.sanger.ac.uk/metadata/control_set_small.all_consequences_with_gene_and_csq.tsv.bgz).

To export VEP consequences of VCF files,
you need to provide the correct VCF header
via the option `stage4 -> annotate_cq_rf -> csq_header`.
To create it, run:

```bash
bcftools view -h dataset.vep-annotated.vcf.gz > dataset..csq_header.txt'
```
The consequence export works only if **both** `general->metadata->vep_consequences`
and `stage4 -> annotate_cq_rf -> csq_header` are specified.

```shell
python stage4_genotype_qc/gqc2_apply_range_of_hard_filters.py
```

#### (Optional) Sex chromosome-specific filtering

By default, the script applies different filtering thresholds for sex chromosomes (X and Y).
For males on chrX non-PAR and chrY regions,
the depth (DP) threshold is halved (since these regions are hemizygous).
This behavior is controlled by the `sex_chromosome_specific_filtering` option
in the `apply_hard_filters` section of the config file.
Set it to `False` to use the same filtering for all chromosomes.

### Export the filtered variants to VCFs

Script `3a` tags all variations with the corresponding filter (relaxed, medium, stringent)
and saves the resulting data to VCF files.
By default, all variants are exported "as is".
You can filter out variants that do not pass relaxed filter by
specifying `drop_failed: False` option for the step 4.3a in the config file.


```shell
python stage4_genotype_qc/gqc3a_export_vcfs_range_of_hard_filters.py
```

Alternatively, to export VCFs with only a single filter level applied,
use the `3b` version of the script.
You can specify which filter level to apply using the `--filter-level` option
(choices: `relaxed`, `medium`, `stringent`; default: `stringent`):

```shell
python stage4_genotype_qc/gqc3b_export_vcfs_single_filter.py --filter-level stringent
```

### Information stored in the VCF data

The resulting VCF contains the following data:
* `FILTER` column tags, describing whether the variant passes stringent, medium, or relaxed filters.
  *Note: The export step removes all variants not passing the relaxed filter.
  Therefore, you won't see any variants with filter tag `relaxed_fail`*
* The full set of variant-level annotations form the original VCFs
* The RF bin and RF score are added to the `INFO` column.
* VEP consequences annotations are added to the `CSQ` entry if were provided on the genotype filtering step
* The full set of genotype-related metrics from original VCfs: DP, GQ, AB
* Genotype-level tags `stringent_filters`, `medium_filters`, and `relaxed_filters`
  describing whether the **genotype** passes the relaxed, medium, or stringent filters.
  *Note: If a variant fails the variant-level filters (like RF bin or call_rate), all genotypes are tagged as `fail`*


### (Optional) - calculate-per-sample statistics and variant statistics
If you want to additionally evaluate the variant QC and genotypeQC statistics,
use the following script:

```shell
python stage4_genotype_qc/gqc5_counts_per_sample.py
```

The script makes the summary table of variations by variation consequence
and saves it the `annotations` by name defined in the `cq_stats_df` config section.
Please note that one variant can be annotated by more than one consequence,
so the sum of stats is not equal to the total number of variations.

Also, the script outputs in the console the variant counts per sample

This step uses the same VEP consequences file specified in `general->metadata->vep_consequences`
If you set it to `null`, the consequence annotation will be skipped,
and per-consequence counts will not be calculated.

The VEP consequences file should be in the following format:

```
chr10   100199947   rs367984062 A   C   intron_variant          intron_variant      GENE1   HGNC:1234
chr10   100199976   rs774723210 G   A   missense_variant        missense_variant    GENE2   HGNC:5678
chr10   100200004   .           C   A   missense_variant        missense_variant    GENE3   HGNC:9012
chr10   100200012   rs144642900 C   T   missense_variant        missense_variant    GENE4   HGNC:3456
chr1    100200019   .           C   A   stop_gained&splice_region_variant  stop_gained  GENE5   HGNC:7890
```


### Plot mutations spectra

After completion of the QC process, run the mutation spectra calculation
and validate that results match the expected distribution.

To do it, run the calculation script:

```shell
python stage4_genotype_qc/gqc4_mutation_spectra_afterqc.py
```

The script saves the plot in the HTML file specified under the
`plot_mutation_spectra_afterqc`:`mut_spectra_path` config section.
