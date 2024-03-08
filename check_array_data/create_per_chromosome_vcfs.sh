SHARD_VCFS="/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/BiB/gatk_vcfs"
VCF_LISTS="/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/BiB/per_chromosome_vcf_lists"
CHROM_VCFS="/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/BiB/gatk_vcfs_per_chromosome"
CMDS="/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/BiB/concat_vcfs.cmds"

mkdir -p $VCF_LISTS $CHROM_VCFS
rm -f $CMDS

for chrom in {1..22}
do
  VCF_LIST="$VCF_LISTS/chr${chrom}.txt"
  find $SHARD_VCFS -name "*chr${chrom}*.gz" | sort -V > $VCF_LIST
  echo bcftools concat -n --file-list $VCF_LIST \| bcftools view -t "chr${chrom}" -Ob -o "${CHROM_VCFS}/chr${chrom}.bcf" >> $CMDS
done

module add common-apps/bcftools/1.19
wr add -f $CMDS -i "ip13.bcftools_concat.bib" --req_grp "bcftools.concat" --cpus 2 --limit_grps "bcftools.concat"

rm $CMDS
# execute after all wr jobs are done
# rm -r $VCF_LISTS
