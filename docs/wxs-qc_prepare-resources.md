# WxS-QC resources howto

This document describes resource files used by WxS-QC pipeline and how to obtain them.
We expect that you have already set up a work machine / cluster and
cloned the pipeline repository to it.

All resources used by WxS-QC pipeline stored in two folders:
`resources` and `training_sets`.

## 1000 Genomes sample data

The SampleQC stage of the WxS-QC pipeline uses data from 1000 genomes
to run clustering and identify superpopulations.

Briefly, you need to do the following:

1. Download per-chromosome 1000G dataset. HGI uses release taken from here:
   ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/
2. (Optional) - run BCFTools to remove structural variations and keep only SNVs and small indels.
3. Put the data in the folder specified under `onekg_resource_dir` in the `general` config section (see below).

### Prerequisites
- Access to a Linux/Unix environment
- `wget` installed
- `bcftools` installed
- `bgzip` and `tabix` installed (you can install it with `samtools` package)
- (optional) Cluster with LSF job scheduling system (`bsub` commands)
- Sufficient storage space for genomic data

### Set up the directory structure

```bash
export KG_DIR=/path/to/your/dataset/folder/resources/1kg_vcf
mkdir -p $KG_DIR
cd $KG_DIR
```

### Download data

The `1kg_download` script downloads VCF files from the 1000 Genomes Project:

```bash
/wes/qc/root/folder/scripts/1kg_download
```

### Remove structural variants
The downloaded 1000 genomes VCFs contain some structural variation.
We need to remove it and keep only SNVs and short indels.

The `1kg_remove_sv` script removes structural variants from VCF files.
You need to run the script from the 1000 genomes data folder with the chromosome name as argument.

```bash
# To process all chromosomes sequentially
for chr in {1..22} X Y; do /wes/qc/root/folder/scripts/1kg_remove_sv $chr; done
```

If you work on a cluster with the installed IBM LSF job scheduling system,
you can use it to submit your jobs.
The script automatically uses LSF job array number for chromosome names

```bash
# Chromosome X
bsub -q small -n 4 /wes/qc/root/folder/scripts/1kg_remove_sv X

# Chromosome Y
bsub -q small -n 4 /wes/qc/root/folder/scripts/1kg_remove_sv Y

# Pick all numeric chromosomes from the job array ID
bsub -q small -n 4 -J [1-22] /wes/qc/root/folder/scripts/1kg_remove_sv
```

### Remove original files

After checking that all files were processed correctly,
you can remove the original 1000 genomes files:

```bash
rm 1kGP_high_coverage_Illumina.chr*.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
```

### (Optional) Reduce data size for exome-only analysis

The manual above prepares the full-genome dataset, containing variations across the whole human genome,
If you work only with exome data, you can significantly reduce the data size,
keeping only exome regions.

First, you need the BED/GFF/GTF file containing coordinates of exome regions.

If you have the BED file for your exome enrichment kit,
we highly recommend you to use it.

If you don't have any, you can download GTF from _Gencode_ and keep only exon regions from it.

First, download all reference data and create genome size file from reference genome.
You need the `samtools` package to generate index containing genome sizes.

```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz

gunzip GRCh38.primary_assembly.genome.fa.gz
samtools faidx GRCh38.primary_assembly.genome.fa
```

Now we extract only exons, add +50 bp to capture flanking regions and merge overlapping regions.
For this work you need the `bedtools` package to manipulate genome regions.

```bash
gunzip -c gencode.v47.annotation.gtf.gz | awk '{FS="\t";OFS="\t"} $3=="exon"' | bedtools slop -g GRCh38.primary_assembly.genome.fa.fai -b 50 | bedtools sort -g GRCh38.primary_assembly.genome.fa.fai | bedtools merge > GRCh38.gencode.v47.exome.bed
```

Finally, extract only exome regions from 1000 genomes VCFs:

```bash
for chr in {1..22} X Y; do tabix 1kGP_chr${chr}.nosv.vcf.gz; done

for chr in {1..22} X Y; do bcftools view -R GRCh38.gencode.v47.exome.bed 1kGP_chr${chr}.nosv.vcf.gz -Oz > 1kGP_chr${chr}.nosv.exome.vcf.gz; done
```

Now you can copy all VCFs with the `nosv.exome.vcf.gz` suffix to a new folder,
and use it on step 1 of the resource preparation stage.

## gnomAD and other resources

The Variant QC part of the pipeline uses population frequencies from the
[gnomAD project](https://gnomad.broadinstitute.org/)
to find _de novo_ variations.
Technically, for this step you can use the original gnomAD exome/genome data.
However, the full-size gnomAD dataset is very big, so we recommend you use
a reduced version, containing only global population frequencies.

Also, the pipeline uses a set of high-quality variations to train the random forest model.

The easiest way to obtain gnomAD frequencies and all other resources
is to download it from the WSI server:

```bash
wget https://wes-qc-data.cog.sanger.ac.uk/WxS-QC_resources.tar
```

You can copy or symlink the
`resources` folder to your data analysis folder (see the main howto).

_The `resources` folder also contains a small subset of 1000-Genomes data.
However, this set is test-only, and for production run
you should download the full-sized 1000-Genomes dataset._


### Using original gnomAD data

If you want to use your own gnomAD data (for example, for genome frequencies),
  you need to manually download it from https://gnomad.broadinstitute.org/downloads
  (use the _Sites Hail Table_ version),
  place the path to the table in the config file section `prepare_gnomad_ht -> input_gnomad_htfile`,
  and run the script to make a reduced version:
  ```shell
  spark-submit 0-resource_preparation/3-prepare-gnomad-table.py
  ```

## Resources description:

### `resources` folder

* `igsr_samples.tsv` -- known super populations for 1000 genomes dataset.
  Available here: https://www.internationalgenome.org/data-portal/sample (press the 'Download the list' button)
* `long_ld_regions.hg38.bed` -- BED file containing long-range linkage disequilibrium regions for the genome version hg38
  The regions were obtained from the file `high-LD-regions-hg38-GRCh38.bed` in **plinkQC** github repo:
  (https://github.com/cran/plinkQC/blob/master/inst/extdata/high-LD-regions-hg38-GRCh38.bed).
  These coordinates are results of `liftOver` transferring original coordinates from the genome version hg36 to hg38.
  Original coordinates are provided in supplementary files of the article
  **Anderson, Carl A., et al. "Data quality control in genetic case-control association studies."
  Nature protocols 5.9 (2010): 1564-1573. DOI: 10.1038/nprot.2010.116**
* `HG001_GRCh38_benchmark.interval.illumina.vcf.gz` with `tabix` index -- High-confidence variations for GIAB HG001 sample.
  Available here: https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
* `HG001_GRCh38_benchmark.all.interval.illumina.vep.info.txt` - VEP annotations for GIAB HG001 sample. Can generate
* `gnomad.exomes.r4.1.freq_only.ht` - reduced version of **gnomAD** data containing only global population frequencies
* `1000G_phase1.snps.high_confidence.hg38`, `1000G_omni2.5.hg38`,
  `hapmap_3.3.hg38`, `Mills_and_1000G_gold_standard.indels.hg38` -
   the set of high-confident variations from the corresponding projects.
   The corresponding VCFs are available from the GATK resource bundle:
   (https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle).
   - Google Cloud: https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/
   - FTP: http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg38/
