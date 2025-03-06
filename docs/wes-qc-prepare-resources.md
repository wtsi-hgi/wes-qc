# How to prepare resource data for WES-QC pipeline

## Preparing the 1000 Genomes sample data

The SampleQC stage of the WES-QC pipline uses data from 1000 genomes
to run clustering and form superpopulations

### Prerequisites
- Access to a Linux/Unix environment
- `wget` installed
- `bcftools` installed
- (optional) LSF job scheduling system (`bsub` commands)
- Sufficient storage space for genomic data


Next these data are used to

To obtain 1000Genomes data, do the following:

1. 1. Set up the directory structure

```bash
export KG_DIR=/path/to/your/dataset/folder/resources/1kg_vcf
mkdir -p $KG_DIR
cd $KG_DIR
```

2. Run the downloading script

The `1kg_download` script downloads VCF files from the 1000 Genomes Project:

```bash
/wes_qc_root_folder/scripts/1kg_download
```

3. Remove structural variants

The `1kg_remove_sv` script processes the downloaded VCF files to remove structural variants.
You need to run the script from the 1kg data folder with the chromosome name

```bash
# To process all chromosomes sequentially
for chr in {1..22} X Y; do
    /wes_qc_root_folder/scripts/1kg_remove_sv $chr
done
```

If you work on a cluster with the installed IBM LSF job scheduling system,
you can use it to submit your jobs as follows.
The script automatically uses LSF job array number for chromosome names

```bash
# Chromosome X
bsub -q small -n 4 /wes/qc/root/folder/scripts/1kg_remove_sv X

# Chromosome Y
bsub -q small -n 4 /wes/qc/root/folder/scripts/1kg_remove_sv Y

# Pick all numeric chromosomes from the job array ID
bsub -q small -n 4 -J [1-22] /wes/qc/root/folder/scripts/1kg_remove_sv
```

4. Remove original files

After checking that all files were processed correctly,
you can remove the original 1000 genomes files:

```bash
rm 1kGP_high_coverage_Illumina.chr*.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
```
