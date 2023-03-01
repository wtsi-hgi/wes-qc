#!/bin/bash

plinkdata="/lustre/scratch126/humgen/projects/birth_cohort_wes/mcs/SNP_array_genotype/full_data_received_201222/MCS_QC/MCS_autosome"
outvcf=/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/mcs/check_array_genotypes/plink_vcf_201222
plink=/software/hgi/installs/plink/dev-05-03-2022/plink2

${plink} --bfile ${plinkdata} --recode vcf bgz id-paste=iid --out ${outvcf}

tabix -p vcf ${outvcf}.vcf.gz
