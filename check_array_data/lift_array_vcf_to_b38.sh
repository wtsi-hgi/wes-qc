#!/bin/bash
set -o pipefail on 

# infile
arrayvcf=/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/plink_vcf.vcf.gz
# wrapper script for lift over
liftoverscript=/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/tools/lift-over
# resources
b37fa=/lustre/scratch118/humgen/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa
b38fa=/lustre/scratch118/humgen/resources/ref/Homo_sapiens/GRCh38_15/Homo_sapiens.GRCh38_15.fa
chain=/lustre/scratch118/humgen/resources/liftover/hg19ToHg38.over.chain
# output
prefix=/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/liftover
outfile=/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/plink_vcf_b38_liftover.vcf.gz

outfile_l=${outfile}_lifted
outfile_header_fix=${outfile}_hf

#liftOver path
liftoverpath=/software/hgi/installs/liftOver

export PERL5LIB=/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/tools/:$PERL5LIB

${liftoverscript} -i ${arrayvcf} -p ${prefix} -o ${b37fa} -n ${b38fa} -c ${chain} -l ${liftoverpath} > ${outfile_l}

sed s'/ID=/ID=chr/' ${oufile_l} > ${outfile_header_fix}
bcftools sort -Oz -o ${outfile} ${outfile_header_fix}

tabix -p vcf ${outfile}

#rm ${outfile_l}
#rm ${outfile_header_fix}