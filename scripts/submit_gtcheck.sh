#!/bin/bash

#BSUB -q long
#BSUB -o /lustre/scratch123/hgi/projects/birth_cohort_wes/qc/BiB/gtcheck/submit_gtcheck.%J.o
#BSUB -e /lustre/scratch123/hgi/projects/birth_cohort_wes/qc/BiB/gtcheck/submit_gtcheck.%J.e
#BSUB -M 20000
#BSUB -R "select[mem>20000] rusage[mem=20000]"

module add common-apps/bcftools/1.19
bcftools gtcheck \
  -S gt:/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/BiB/gtcheck/microarray_data.samples.txt \
  -S qry:/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/BiB/gatk_vcfs_sites_in_array/all_chroms_isec.samples.txt \
  --n-matches 50 \
  -g /lustre/scratch123/hgi/projects/birth_cohort_wes/qc/BiB/gtcheck/microarray_data.b38_liftover.bcf \
  /lustre/scratch123/hgi/projects/birth_cohort_wes/qc/BiB/gatk_vcfs_sites_in_array/all_chroms_isec.bcf \
  > /lustre/scratch123/hgi/projects/birth_cohort_wes/qc/BiB/gtcheck/gtcheck_top_50_matches.txt
