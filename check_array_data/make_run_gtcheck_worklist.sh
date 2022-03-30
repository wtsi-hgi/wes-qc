#!/bin/bash
 
vcflist="/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/gatk_vcf_list.txt"

bcftools="/nfs/users/nfs_p/pd3/bcftools/dist/bin/bcftools"
samplesg="/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/plink_samples.txt"
samplesv="/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/gatk_samples.txt"
plinkvcf="/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/plink_vcf_b38_liftover.bcf"
outdir="/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/run_gtcheck_all_pairs/"
outfile=${outdir}"gtcheck_out_"${LSB_JOBINDEX}

worklist="/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/run_gtcheck_worklist_all_pairs.txt"

for f in `cat ${vcflist}`
    do
    vcf="/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/gatk_vcfs/"${f}
    pos=${f::-12}
    cmd=${bcftools}" gtcheck -S gt:"${samplesg}" -S qry:"${samplesv}" -g "${plinkvcf}" "${vcf}" > "${outfile}${pos}
    echo '"'${cmd}'"'
    done > ${worklist}