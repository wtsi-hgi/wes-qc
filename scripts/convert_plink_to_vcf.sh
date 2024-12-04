#!/bin/bash

plinkdata=/lustre/scratch123/hgi/mdt1/projects/autozyg/BiB/files_from_Mark/BIBmerged_6
outvcf=/lustre/scratch123/hgi/mdt2/projects/birth_cohort_wes/qc/BiB/gtcheck/microarray_data
plink=/software/hgi/installs/plink/dev-05-03-2022/plink2

${plink} --bfile ${plinkdata} --recode vcf bgz id-paste=fid --out ${outvcf}

tabix -p vcf ${outvcf}.vcf.gz
