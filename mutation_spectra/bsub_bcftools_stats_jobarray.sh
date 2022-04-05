#!/bin/bash
#run bcftools ststs as a job array - for each job index run on all 24 chromsomes sequentially as each job is failry short

sample_list=/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/gatk_samples.txt
wdir=$(pwd)

script=${wdir}/run_bcftools_stats.sh
echo ${script}
exit 0

sample=$(sed -n -e ${LSB_JOBINDEX}p ${sample_list})

for chr in {{1..22},X,Y}
do
script ${sample} ${chr}
done

