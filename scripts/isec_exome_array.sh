ARRAY_VCF="/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/BiB/gtcheck/microarray_data.b38_liftover.bcf"
EXOME_VCFS="/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/BiB/gatk_vcfs_per_chromosome"
ISEC_VCFS="/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/BiB/gatk_vcfs_sites_in_array"
CMDS="/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/BiB/isec_vcfs.cmds"

rm -f "${CMDS}"
mkdir -p "${ISEC_VCFS}"

for chrom in {1..22}
do
  EXOME_VCF="${EXOME_VCFS}/chr${chrom}.bcf"
  echo bcftools isec -n =2 -w 1 -Ob -o "${ISEC_VCFS}/chr${chrom}_isec.bcf" $EXOME_VCF $ARRAY_VCF >> $CMDS
done

module add common-apps/bcftools/1.19
wr add -f $CMDS -i "ip13.bcftools_isec.bib" --req_grp "bcftools.isec" --limit_grps "bcftools.isec"

rm $CMDS

bcftools concat -f to_concat.txt -n -Ob -o /lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/gatk_vcfs_sites_in_plink/all_chroms_isec.bcf.gz
