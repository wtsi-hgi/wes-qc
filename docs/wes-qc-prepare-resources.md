# How to prepare resource data for WES-QC pipeline

## Preparing the 1000 Genomes sample data

The SampleQC stage of the WES-QC pipline uses data from 1000 genomes
to run clustering and form super populations.

### Prerequisites
- Access to a Linux/Unix environment
- `wget` installed
- `bcftools` installed
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
/wes_qc_root_folder/scripts/1kg_download
```

### Remove structural variants
The downloaded 1000 genomes VCFs contain some structural variation.
We nned to remove it and keep only SNVs and short indels.

The `1kg_remove_sv` script removes structural variants from VCF files.
You need to run the script from the 1000 genomes data folder with the chromosome name as argument.

```bash
# To process all chromosomes sequentially
for chr in {1..22} X Y; do /wes_qc_root_folder/scripts/1kg_remove_sv $chr; done
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

The manual above prepares the full-genome dataset, containing variations across all genome,
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

Finally, we extract only exome regions from 1000 genomes VCFs:

```bash
for chr in {1..22} X Y; do tabix 1kGP_chr${chr}.nosv.vcf.gz; done

for chr in {1..22} X Y; do bcftools view -R GRCh38.gencode.v47.exome.bed 1kGP_chr${chr}.nosv.vcf.gz -Oz > 1kGP_chr${chr}.nosv.exome.vcf.gz; done
```

Now you can copy all VCFs with the `nosv.exome.vcf.gz` suffix to a new folder,
and use it on step 1 of the resource preparation stage.
