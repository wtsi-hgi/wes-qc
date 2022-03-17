#!/bin/bash

plinkdata=/lustre/scratch119/humgen/projects/birth_cohort_wes/alspac/existing_data/ARRAY_2020-09-03/data
outvcf=/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/plink_vcf
plink=/software/hgi/installs/plink/dev-05-03-2022/plink2

${plink} --bfile ${plinkdata} --recode vcf --out ${outvcf}

bgzip -c ${outvcf}.vcf > ${outvcf}.vcf.gz

tabix -p vcf ${outvcf}.vcf.gz
