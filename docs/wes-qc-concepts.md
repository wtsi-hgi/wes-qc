---
tags:
  - Concept
  - Compute
  - Farm
  - Pipelines
---

# WES QC

The aim of WES QC is to produce a clean dataset removing problem samples, and excluding variants and genotypes that are liekely to be false positives. Each dataset is different and the QC requires careful tailoring and evaluation at each step to ensure that it is appropriate for the dataset being worked on.

## Overview

WES QC consists of three main steps: Sample QC, variant QC and genotype QC. It is written using (Hail)[https://hail.is/] and much of it is based on methods used by gnomAD:

https://gnomad.broadinstitute.org/news/2018-10-gnomad-v2-1/  
https://gnomad.broadinstitute.org/news/2019-10-gnomad-v3-0/  
https://gnomad.broadinstitute.org/news/2020-10-gnomad-v3-1-new-content-methods-annotations-and-data-availability/


## Sample QC

The aim of sample QC is to remove poor quality samples. These may be samples where the sequencing is poor quality, these may be contaminated samples or these may be samples for which the identity is in question.

### Stratified sample QC

The superpopulation each sample belongs to is inferred using PCA. Outliers are identified within each population for each metric tested: 
- Number of SNPs  
- Number of deletions  
- Number of insertions  
- Insertion/deletion ratio  
- Heterozygosity rate  
- Heterozygous/homozygous ratio  
- Number of transitions  
- Number of transversions  
- Transition/transversion ratio.

### Contaminated samples

freemix

### Sample identity

gtcheck

## Variant QC

## Genotype QC
