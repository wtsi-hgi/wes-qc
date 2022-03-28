#!/bin/bash
 
# cmd=$(sed -n -e ${LSB_JOBINDEX}p /lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/run_gtcheck_worklist.txt)

# echo $cmd
# $cmd

for c in `cat /lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/run_gtcheck_worklist.txt`
    do
    bsub -qnormal -o /lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/farmout/gtcheck_array.o.%J -e /lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/farmout/gtcheck_array.e.%J -R"select[mem>10000] rusage[mem=10000]" -M10000 $c
    done