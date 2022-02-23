chr=$LSB_JOBINDEX
if [ $chr = 23 ]
then
    chr='X'
elif [ $chr = 24 ]
then
    chr='Y'    
fi

kgdir="/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/1kg_pca_hail/1kg_data/"
outdir="/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/resources/1kg_vcfs_filtered_by_wes_baits/"

regions="/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/resources/targets.txt"

vcfin=${kgdir}"CCDG_13607_B01_GRM_WGS_2019-02-19_chr"${chr}".recalibrated_variants.vcf.gz"
vcfout=${outdir}"1kg_wes_regions_chr"${chr}".vcf.gz"

bcftools view -R ${regions} -Oz -o ${vcfout} ${vcfin} 
bcftools index ${vcfout}