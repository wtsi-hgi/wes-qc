#!/bin/bash

# infile
arrayvcf=/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/plink_vcf.vcf.gz
# wrapper script for lift over
liftoverscript=/nfs/users/nfs_p/pd3/cvs/vr-codebase/scripts/lift-over
# resources
b37fa=/lustre/scratch119/humgen/projects/ddd/resources/v1.2/hs37d5.fasta
b38fa=/lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/hs38DH.fa
chain=/lustre/scratch118/humgen/resources/liftover/hg19ToHg38.over.chain
# output
prefix=/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/liftover
outfile=/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/plink_vcf_b38_liftover.vcf.gz

${liftoverscript} -i ${arrayvcf} -p ${prefix} -o ${b37fa} -n ${b38fa} -c ${chain} | bgzip -c > ${outfile}

tabix -p vcf ${outfile}