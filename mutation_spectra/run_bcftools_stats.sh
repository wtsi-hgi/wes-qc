#!/bin/bash

sample=$1
chr=$2

bcftools=/nfs/users/nfs_p/pd3/bcftools/dist/bin/bcftools
bcfdir=/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/gatk_vcfs/per_chromosome/
outdir=/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/mutation_spectra/

${bcftools} view -s ${sample}  --trim-alt-alleles -Ou ${bcfdir}chr${chr}.bcf | ${bcftools} stats  -i 'N_ALT>1' > ${outdir}${sample}_chr${chr}.stats