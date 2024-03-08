#!/bin/bash
#run bcftools stats as a job array - for each job index run on all 24 chromosomes sequentially as each job is fairly short

wdir=$(pwd)
sample_list=/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/BiB/gatk_vcfs_sites_in_array/all_chroms_isec.samples.txt
cmds_file=${wdir}/mutation_spectra.cmds

script=${wdir}/run_bcftools_stats.sh
echo ${script}

parallel --dry-run ${script} {} :::: $sample_list > ${cmds_file}

#wr add -f ${cmds_file} -i "ip13.bcftools_stats.bib" -g "bcftools.stats" -l "bcftools.stats:500"

#rm ${cmds_file}
