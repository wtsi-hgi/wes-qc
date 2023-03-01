#!/bin/bash
set -o pipefail on 

# infile
arrayvcf=/lustre/scratch123/hgi/mdt2/projects/birth_cohort_wes/qc/mcs/check_array_genotypes/plink_vcf_201222.vcf.gz
# wrapper script for lift over
liftoverscript=/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/tools/lift-over
# resources
b37fa=/lustre/scratch125/humgen/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa
b38fa=/lustre/scratch125/humgen/resources/ref/Homo_sapiens/GRCh38_15/Homo_sapiens.GRCh38_15.fa
chain=/lustre/scratch125/humgen/resources/liftover/hg19ToHg38.over.chain
# output
prefix=/lustre/scratch123/hgi/mdt2/projects/birth_cohort_wes/qc/mcs/check_array_genotypes/liftover
outfile=/lustre/scratch123/hgi/mdt2/projects/birth_cohort_wes/qc/mcs/check_array_genotypes/plink_vcf_201222_b38_liftover.bcf

outfile_l=${outfile}_lifted
outfile_header_fix=${outfile}_hf

#liftOver path
liftoverpath=/software/hgi/installs/liftOver

export PERL5LIB=/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/tools/:$PERL5LIB

${liftoverscript} -i ${arrayvcf} -p ${prefix} -o ${b37fa} -n ${b38fa} -c ${chain} -l ${liftoverpath} | bgzip > ${outfile_l}

bcftools reheader --fai ${b38fa}.fai ${outfile_l} > ${outfile_header_fix}
bcftools sort -Ob -o ${outfile} ${outfile_header_fix}
bcftools index ${outfile}

#rm ${outfile_l}
#rm ${outfile_header_fix}

bsub -qlong -o farmout/gtcheck.o.%J -e farmout/gtcheck.e.%J -R"select[mem>20000] rusage[mem=20000]" -M20000 "/nfs/users/nfs_p/pd3/bcftools/dist/bin/bcftools gtcheck -S gt:/lustre/scratch123/hgi/mdt2/projects/birth_cohort_wes/qc/mcs/check_array_genotypes/plink_vcf_201222.samples -S qry:/lustre/scratch123/hgi/mdt2/projects/birth_cohort_wes/qc/mcs/check_array_genotypes/gatk_vcfs_sites_in_plink/all_chroms_isec.samples --n-matches 50 -g /lustre/scratch123/hgi/mdt2/projects/birth_cohort_wes/qc/mcs/check_array_genotypes/plink_vcf_201222_b38_liftover.bcf /lustre/scratch123/hgi/mdt2/projects/birth_cohort_wes/qc/mcs/check_array_genotypes/gatk_vcfs_sites_in_plink/all_chroms_isec.bcf > /lustre/scratch123/hgi/mdt2/projects/birth_cohort_wes/qc/mcs/check_array_genotypes/gtcheck_top_50_matches.txt"