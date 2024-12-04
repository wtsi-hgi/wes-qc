#!/bin/bash

cmd=$(sed -n -e ${LSB_JOBINDEX}p /lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/run_gtcheck_worklist.txt)

echo $cmd
$cmd
