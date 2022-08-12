#samples sequenced and variant caLLed at Sanger and at Broad
#50 of the samples sequenced and called at both sites have a 3-fold excess of SNPs and were sequenced twice at Sanger
#for these samples compare broad, original sanger data and resequenced sanger data

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

    #load mts
    broad_mtfile = mtdir + "gatk_calls_from_broad.mt"
    sanger_resequenced_mt_file = mtdir + 'resequenced_samples_gatk_unprocessed.mt'
    sanger_high_snp_samples_mtfile = mtdir + "gatk_unprocessed_high_snp_samples.mt"
    broad_mt = hl.read_matrix_table(broad_mtfile)
    sanger_reseq_mt = hl.read_matrix_table(sanger_resequenced_mt_file)
    sanger_orig_mt = hl.read_matrix_table(sanger_high_snp_samples_mtfile)

    #filter to just the samples we want
    broad_id_file = "file:///lustre/scratch123/qc/compare_broad_sanger_vcfs/broad_ids_high_snp_in_both_cohorts.txt"
    sanger_id_file = "file:///lustre/scratch123/qc/compare_broad_sanger_vcfs/sanger_ids_high_snp_in_both_cohorts.txt"
    sanger_sample_ht = hl.import_table(sanger_id_file, delimiter="\t").key_by('s')
    broad_sample_ht = hl.import_table(broad_id_file, delimiter="\t").key_by('s')
    sanger_reseq_mt = sanger_reseq_mt.filter_cols(hl.is_defined(sanger_sample_ht[sanger_reseq_mt.s]))
    sanger_orig_mt = sanger_orig_mt.filter_cols(hl.is_defined(sanger_sample_ht[sanger_orig_mt.s]))
    broad_mt = broad_mt.filter_cols(hl.is_defined(broad_sample_ht[broad_mt.s]))
    
    #split multi
    broad_mt_split = hl.split_multi_hts(broad_mt)
    sanger_reseq_mt_split = hl.split_multi_hts(sanger_reseq_mt)
    sanger_orig_mt_split = hl.split_multi_hts(sanger_orig_mt)

    #save split multi and reload
    broad_split_multi_mtfile = mtdir + "gatk_calls_from_broad_split_multi_resequenced_samples.mt"
    broad_mt_split.write(broad_split_multi_mtfile, overwrite = True)
    broad_mt = hl.read_matrix_table(broad_split_multi_mtfile)
    sanger_reseq_split_multi_mtfile = mtdir + "resequenced_samples_gatk_unprocessed_split_multi.mt"
    sanger_reseq_mt_split.write(sanger_reseq_split_multi_mtfile, overwrite = True)
    sanger_reseq_mt = hl.read_matrix_table(sanger_reseq_split_multi_mtfile)
    sanger_orig_split_multi_mtfile = mtdir + "gatk_unprocessed_high_snp_samples_split_multi.mt"
    sanger_orig_mt_split.write(sanger_orig_split_multi_mtfile, overwrite = True)
    sanger_orig_mt = hl.read_matrix_table(sanger_orig_split_multi_mtfile)

    #remove any that are only ref alleles
    sanger_reseq_mt = hl.variant_qc(sanger_reseq_mt)
    sanger_orig_mt = hl.variant_qc(sanger_orig_mt)
    broad_mt = hl.variant_qc(broad_mt)
    sanger_reseq_mt = sanger_reseq_mt.filter_rows(sanger_reseq_mt.variant_qc.n_non_ref == 0, keep = False)
    sanger_orig_mt = sanger_orig_mt.filter_rows(sanger_orig_mt.variant_qc.n_non_ref == 0, keep = False)
    broad_mt = broad_mt.filter_rows(broad_mt.variant_qc.n_non_ref == 0, keep = False)

    #extract variants in sanger but not broad and vice verse
    #convert sanger and broad mt rows into ht with select_rows (will this be ht or mt?)
    broad_vars = broad_mt.rows()
    sanger_reseq_vars = sanger_reseq_mt.rows()
    sanger_orig_vars = sanger_orig_mt.rows()
    #select rows in sanger which are not in broad and vice versa
    sanger_orig_not_broad_mt = sanger_orig_mt.anti_join_rows(broad_vars)
    sanger_reseq_not_broad_mt = sanger_reseq_mt.anti_join_rows(broad_vars)
    broad_not_sanger_reseq_mt = broad_mt.anti_join_rows(sanger_reseq_vars)
    broad_not_sanger_orig_mt = broad_mt.anti_join_rows(sanger_orig_vars)

    #save sanger and broad mts and mts unique to each set
    broad_mt.write(broad_split_multi_mtfile, overwrite = True)
    sanger_reseq_mt.write(sanger_reseq_split_multi_mtfile, overwrite = True)
    sanger_orig_mt.write(sanger_orig_split_multi_mtfile, overwrite = True)
    sanger_reseq_only_mtfile = mtdir + "sanger_reseq_variants_not_in_broad.mt"
    sanger_reseq_not_broad_mt.write(sanger_reseq_only_mtfile, overwrite = True)
    sanger_orig_only_mtfile = mtdir + "sanger_orig_variants_not_in_broad.mt"
    sanger_orig_not_broad_mt.write(sanger_orig_only_mtfile, overwrite = True)
    broad_not_sanger_reseq_mtfile = mtdir + "broad_variants_not_in_sanger_reseq.mt"
    broad_not_sanger_reseq_mt.write(broad_not_sanger_reseq_mtfile, overwrite = True)
    broad_not_sanger_orig_mtfile = mtdir + "broad_variants_not_in_sanger_orig.mt"
    broad_not_sanger_orig_mt.write(broad_not_sanger_orig_mtfile, overwrite = True)

    #now repeat for after variant QC
    #orig = runhash b1704ad6, cut off snvs = 17, indels = 57
    #reseq = runhash 746dff9a, cut off SNVs = 89, indels = 84
    rf_htfile_orig = rf_dir +  "/b1704ad6/_gnomad_score_binning_tmp.ht"
    rf_ht_orig = hl.read_table(rf_htfile_orig)
    sanger_orig_mt = sanger_orig_mt.annotate_rows(
        info=sanger_orig_mt.info.annotate(
        rf_bin=rf_ht_orig[sanger_orig_mt.row_key].bin)
    )
    sanger_orig_mt_rf_pass = sanger_orig_mt.filter_rows( ( (hl.is_snp(sanger_orig_mt.alleles[0], sanger_orig_mt.alleles[1]) & (sanger_orig_mt.info.rf_bin <= 17) )
                                             | (( hl.is_indel(sanger_orig_mt.alleles[0], sanger_orig_mt.alleles[1]) ) &  (sanger_orig_mt.info.rf_bin <= 57)) ) )


    sanger_orig_vars_rf_pass = sanger_orig_mt_rf_pass.rows()
    sanger_orig_only_mt_post_rf = sanger_orig_mt_rf_pass.anti_join_rows(broad_vars)
    broad_not_sanger_orig_post_rf = broad_mt.anti_join_rows(sanger_orig_vars_rf_pass)

    sanger_orig_filtered_mtfile_rf_pass = mtdir + "gatk_calls_from_sanger_orig_samples_in_broad_after_rf.mt"
    sanger_orig_only_mt_post_rf.write(sanger_orig_filtered_mtfile_rf_pass, overwrite = True)
    sanger_orig_only_mtfile_rf_pass = mtdir + "sanger_orig_variants_not_in_broad_after_rf.mt"
    sanger_orig_only_mt_post_rf.write(sanger_orig_only_mtfile_rf_pass, overwrite = True)
    broad_not_sanger_orig_mtfile_rf_pass = mtdir + "broad_variants_not_in_sanger_orig_after_rf.mt"
    broad_not_sanger_orig_post_rf.write(broad_not_sanger_orig_mtfile_rf_pass, overwrite = True)


    rf_htfile_reseq = rf_dir +  "/746dff9a/_gnomad_score_binning_tmp.ht"
    rf_ht_reseq = hl.read_table(rf_htfile_reseq)
    sanger_reseq_mt = sanger_reseq_mt.annotate_rows(
        info=sanger_reseq_mt.info.annotate(
        rf_bin=rf_ht_reseq[sanger_reseq_mt.row_key].bin)
    )
    sanger_reseq_mt_rf_pass = sanger_reseq_mt.filter_rows( ( (hl.is_snp(sanger_reseq_mt.alleles[0], sanger_reseq_mt.alleles[1]) & (sanger_reseq_mt.info.rf_bin <= 89) )
                                             | (( hl.is_indel(sanger_reseq_mt.alleles[0], sanger_reseq_mt.alleles[1]) ) &  (sanger_reseq_mt.info.rf_bin <= 84)) ) )

    sanger_reseq_vars_rf_pass = sanger_reseq_mt_rf_pass.rows()
    sanger_reseq_only_mt_post_rf = sanger_reseq_mt_rf_pass.anti_join_rows(broad_vars)
    broad_not_sanger_reseq_post_rf = broad_mt.anti_join_rows(sanger_reseq_vars_rf_pass)

    sanger_reseq_filtered_mtfile_rf_pass = mtdir + "gatk_calls_from_sanger_reseq_samples_in_broad_after_rf.mt"
    sanger_reseq_only_mt_post_rf.write(sanger_reseq_filtered_mtfile_rf_pass, overwrite = True)
    sanger_reseq_only_mtfile_rf_pass = mtdir + "sanger_reseq_variants_not_in_broad_after_rf.mt"
    sanger_reseq_only_mt_post_rf.write(sanger_reseq_only_mtfile_rf_pass, overwrite = True)
    broad_not_sanger_reseq_mtfile_rf_pass = mtdir + "broad_variants_not_in_sanger_reseq_after_rf.mt"
    broad_not_sanger_reseq_post_rf.write(broad_not_sanger_reseq_mtfile_rf_pass, overwrite = True)


if __name__ == '__main__':
    main()