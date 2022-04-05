#!/bin/bash

sample=$1
chr=$2

bcftools=/nfs/users/nfs_p/pd3/bcftools/dist/bin/bcftools
bcfdir=/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/gatk_vcfs/per_chromosome/
outdir=/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/mutation_spectra/
tmpdir=/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/mutation_spectra/tmp/
tmpvcf=tmpdir${sample}_chr${chr}.bcf

echo "Running with sample:"$1" chromosome:chr"$2

${bcftools} view -s ${sample} --trim-alt-alleles -Ob -o ${tmpbcf} ${bcfdir}chr${chr}.bcf 
${bcftools} index ${tmpvcf}

${bcftools} stats  -i 'N_ALT>1' ${tmpbcf} > ${outdir}${sample}_chr${chr}.stats

