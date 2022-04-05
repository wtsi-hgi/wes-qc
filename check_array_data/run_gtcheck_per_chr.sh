#!/bin/bash

c=$LSB_JOBINDEX

bcftools=/nfs/users/nfs_p/pd3/bcftools/dist/bin/bcftools
samplesg=/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/plink_samples.txt
samplesv=/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/gatk_samples.txt
plinkvcf=/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/plink_vcf_b38_liftover.bcf
outdir=/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/per_chromosome_output/
outfile=${outdir}chr${c}_gtcheck.txt
vcfdir=/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/gatk_vcfs/per_chromosome/
vcf=${vcfdir}chr${c}.bcf

${bcftools} gtcheck -S gt:${samplesg} -S qry:${samplesv} -g ${plinkvcf} ${vcf} > ${outfile}
