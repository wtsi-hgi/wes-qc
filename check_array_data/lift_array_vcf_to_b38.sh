#!/bin/bash

# infile
arrayvcf=/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/plink_vcf.vcf.gz
# wrapper script for lift over
liftoverscript=/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/tools/lift-over
# resources
b37fa=/lustre/scratch118/humgen/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa
#b38fa=/lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/hs38DH.fa
b38fa=/lustre/scratch118/humgen/resources/ref/Homo_sapiens/GRCh38_15/Homo_sapiens.GRCh38_15.fa
chain=/lustre/scratch118/humgen/resources/liftover/hg19ToHg38.over.chain
#chain=/lustre/scratch118/humgen/resources/liftover/grch37tohg38.chain
# output
prefix=/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/liftover
outfile=/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/plink_vcf_b38_liftover.vcf.gz

#liftOver path
liftoverpath=/software/hgi/installs/liftOver

export PERL5LIB=/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/tools/:$PERL5LIB

${liftoverscript} -i ${arrayvcf} -p ${prefix} -o ${b37fa} -n ${b38fa} -c ${chain} -l ${liftoverpath} | bgzip -c > ${outfile}

#tabix -p vcf ${outfile}