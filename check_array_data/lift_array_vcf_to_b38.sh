#!/bin/bash

#BSUB -q normal
#BSUB -o /lustre/scratch123/hgi/mdt2/projects/birth_cohort_wes/qc/BiB/gtcheck/lift_array_vcf_to_b38.%J.o
#BSUB -e /lustre/scratch123/hgi/mdt2/projects/birth_cohort_wes/qc/BiB/gtcheck/lift_array_vcf_to_b38.%J.e
#BSUB -R "select[mem>20000] rusage[mem=20000]"
#BSUB -M 20000
#BSUB -J liftover_BiB

set -e -o pipefail on

module add common-apps/bcftools/1.19
module add common-apps/samtools/1.19

# infile
arrayvcf=/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/BiB/gtcheck/microarray_data.vcf.gz
# wrapper script for lift over
liftoverscript=/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/tools/lift-over
# resources
b37fa=/lustre/scratch125/humgen/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa
b38fa=/lustre/scratch125/humgen/resources/ref/Homo_sapiens/GRCh38_15/Homo_sapiens.GRCh38_15.fa
chain=/lustre/scratch125/humgen/resources/liftover/hg19ToHg38.over.chain
# output
prefix=/lustre/scratch123/hgi/mdt2/projects/birth_cohort_wes/qc/BiB/gtcheck/liftover
outfile=/lustre/scratch123/hgi/mdt2/projects/birth_cohort_wes/qc/BiB/gtcheck/microarray_data.b38_liftover.bcf

outfile_l=${outfile}_lifted
outfile_header_fix=${outfile}_hf

#liftOver path
liftoverpath=/software/hgi/installs/liftOver

export PERL5LIB=/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/tools/:$PERL5LIB

${liftoverscript} -i ${arrayvcf} -p ${prefix} -o ${b37fa} -n ${b38fa} -c ${chain} -l ${liftoverpath} | bgzip > ${outfile_l}

bcftools reheader --fai ${b38fa}.fai ${outfile_l} -o ${outfile_header_fix}
bcftools sort -Ob -o ${outfile} ${outfile_header_fix}
bcftools index ${outfile}

#rm ${outfile_l}
#rm ${outfile_header_fix}