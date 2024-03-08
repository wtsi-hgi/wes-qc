#!/bin/bash

sample=$1

module add common-apps/bcftools
bcfdir=/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/BiB/gatk_vcfs_per_chromosome
outdir=/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/BiB/mutation_spectra

for chr in {1..22} X Y
do

  outfile=${outdir}/${sample}_${chr}.stats
  bcf=${bcfdir}/chr${chr}.bcf

  echo "Running with sample: $1 chromosome:chr$2"

  bcftools view -s ${sample} --trim-alt-alleles -Ou ${bcf} | bcftools stats -i 'N_ALT>0' > ${outfile}
  bgzip ${outfile}

done
