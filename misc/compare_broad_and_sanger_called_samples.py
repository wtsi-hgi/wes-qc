#samples sequenced and variant caLLed at Sanger and at Broad
#remove samples with a ~3-fold excess of SNPs in the Sanger data and samples where there is any doubt over identify
import hail as hl
import pyspark
from wes_qc.utils.utils import parse_config


def main():
    # set up
    inputs = parse_config()
    mtdir = inputs['matrixtables_lustre_dir']
    rf_dir = inputs['var_qc_rf_dir']

    # initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    broad_mtfile = mtdir + "gatk_calls_from_broad.mt"
    broad_vcf = "file:///lustre/scratch123/qc/compare_broad_sanger_vcfs/broad_vcf_samples_in_sanger.vcf.gz"

    #load broad vcfs into mtfile and save mtfile
    broad_mt = hl.import_vcf(broad_vcf, force_bgz = True)
    broad_mt.write(broad_mtfile, overwrite = True)
    #open sanger mtfile
    samger_mtfile = mtdir + 'gatk_unprocessed.mt'
    sanger_mt = hl.read_matrix_table(samger_mtfile)
    #filter both to only the samples we want
    sanger_samples_file = "file:///lustre/scratch123/qc/compare_broad_sanger_vcfs/sanger_accs_to_analyse_s.txt"
    broad_samples_file = "file:///lustre/scratch123/qc/compare_broad_sanger_vcfs/broad_accs_to_analyse_s.txt"
    sanger_sample_ht = hl.import_table(sanger_samples_file, delimiter="\t").key_by('s')
    broad_sample_ht = hl.import_table(broad_samples_file, delimiter="\t").key_by('s')
    sanger_mt = sanger_mt.filter_cols(hl.is_defined(sanger_sample_ht[sanger_mt.s]))
    broad_mt = broad_mt.filter_cols(hl.is_defined(broad_sample_ht[broad_mt.s]))
    #split multi
    sanger_mt = hl.split_multi_hts(sanger_mt)
    broad_mt = hl.split_multi_hts(broad_mt)
    #save split multi and reload
    broad_split_multi_mtfile = mtdir + "gatk_calls_from_broad_split_multi.mt"
    sanger_split_multi_mtfile = mtdir + "gatk_unprocessed_split_multi.mt"
    broad_mt.write(broad_split_multi_mtfile)
    sanger_mt.write(sanger_split_multi_mtfile)
    sanger_mt = hl.read_matrix_table(sanger_split_multi_mtfile)
    broad_mt = hl.read_matrix_table(broad_split_multi_mtfile)
    #remove any that are only ref alleles
    sanger_mt = hl.variant_qc(sanger_mt)
    broad_mt = hl.variant_qc(broad_mt)
    sanger_mt = sanger_mt.filter_rows(sanger_mt.variant_qc.n_non_ref == 0, keep = False)
    broad_mt = broad_mt.filter_rows(broad_mt.variant_qc.n_non_ref == 0, keep = False)
    #extract variants in sanger but not broad
    #convert sanger and broad mt rows into ht with select_rows (will this be ht or mt?)
    broad_vars = broad_mt.rows()
    sanger_vars = sanger_mt.rows()
    #select rows in sanger which are not in broad (can do this with select those in sanger that are in broad and keep=false) try anti_join_rows https://hail.is/docs/0.2/hail.MatrixTable.html#hail.MatrixTable.anti_join_rows
    sanger_only_mt = sanger_mt.anti_join_rows(broad_vars)
    broad_only_mt = broad_mt.anti_join_rows(sanger_vars)
    #save sanger and broad mts and unique to sanger ht
    broad_filtered_mtfile = mtdir + "gatk_calls_from_broad_samples_in_sanger.mt"
    broad_mt.write(broad_filtered_mtfile, overwrite = True)
    sanger_filtered_mtfile = mtdir + "gatk_calls_from_sanger_samples_in_broad.mt"
    sanger_mt.write(sanger_filtered_mtfile, overwrite = True)
    sanger_only_mtfile = mtdir + "sanger_variants_not_in_broad.mt"
    sanger_only_mt.write(sanger_only_mtfile, overwrite = True)
    broad_only_mtfile = mtdir + "broad_variants_not_in_sanger.mt"
    broad_only_mt.write(broad_only_mtfile, overwrite = True)
    
    #now repeat for after variant QC
    rf_htfile = rf_dir +  "/1ed4bbbc/_gnomad_score_binning_tmp.ht"
    rf_ht = hl.read_table(rf_htfile)
    sanger_mt = sanger_mt.annotate_rows(
        info=sanger_mt.info.annotate(
        rf_bin=rf_ht[sanger_mt.row_key].bin)
    )
    #filter using a pass of bin 37 for SNPs and bin 80 for indels
    sanger_mt_rf_pass = sanger_mt.filter_rows( ( (hl.is_snp(sanger_mt.alleles[0], sanger_mt.alleles[1]) & (sanger_mt.info.rf_bin <= 37) )
                                             | (( hl.is_indel(sanger_mt.alleles[0], sanger_mt.alleles[1]) ) &  (sanger_mt.info.rf_bin <= 80)) ) )


    sanger_vars_rf_pass = sanger_mt_rf_pass.rows()
    sanger_only_mt_post_rf = sanger_mt_rf_pass.anti_join_rows(broad_vars)
    broad_only_mt_post_rf = broad_mt.anti_join_rows(sanger_vars_rf_pass)

    sanger_filtered_mtfile_rf_pass = mtdir + "gatk_calls_from_sanger_samples_in_broad_after_rf.mt"
    sanger_mt_rf_pass.write(sanger_filtered_mtfile_rf_pass, overwrite = True)
    sanger_only_mtfile_rf_pass = mtdir + "sanger_variants_not_in_broad_after_rf.mt"
    sanger_only_mt_post_rf.write(sanger_only_mtfile_rf_pass, overwrite = True)
    broad_only_mtfile_rf_pass = mtdir + "broad_variants_not_in_sanger_after_rf.mt"
    broad_only_mt_post_rf.write(broad_only_mtfile_rf_pass, overwrite = True)
    

if __name__ == '__main__':
    main()