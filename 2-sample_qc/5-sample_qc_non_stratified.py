# perform hail sample QC on raw matrixtable and annotate with superpops and fail count from stratified qc
import resource
import hail as hl
import pyspark
from wes_qc.utils.utils import parse_config


def annotate_mt(raw_mt_file, pop_ht_file, runid_file, annotated_mt_file):
    '''
    Annotate mt with superpopulation and sequencing runid
    :param str raw_mt_file: raw mt file
    :param str pop_ht_file: file with popuklation annotations
    :param str runid_file: metadata file to annotate run ids
    :param str annotated_mt_file: annotated mt file
    '''
    mt = hl.read_matrix_table(raw_mt_file)
    runida_ht = hl.import_table(runid_file, delimiter="\t").key_by('ega')
    mt = mt.annotate_cols(batch=runida_ht[mt.s]['runid'])
    seq_expr = (hl.case()
                .when(mt.s.startswith('EGAN'), 'Sanger')
                .when(mt.s.startswith('Z'), 'Bristol')
                .default("")
                )
    mt = mt.annotate_cols(sequencing_location = seq_expr).key_cols_by('s')
    pop_ht = hl.read_table(pop_ht_file)
    mt = mt.annotate_cols(assigned_pop=pop_ht[mt.s].pop)
    mt.write(annotated_mt_file, overwrite=True)


def run_sample_qc(annotated_mt_file, mt_sample_qc_file, ht_sample_qc_file, annotdir, min_dp: int = 0, min_gq: int = 0, min_vaf: float = 0.0):
    '''
    Run sample QC on raw variants
    :param str annotated_mt_file: annotated mt file 
    :param str mt_sample_qc_file: output MT
    :param str ht_sample_qc_file: output HT, cols only
    :param str annotdir: annotations dir
    :param int min_dp: minimum DP
    :param int min_gq: minimum GQ
    :param flaot min_vaf: minimum VAF
    '''
    mt = hl.read_matrix_table(annotated_mt_file)
    #filter to autosomes only
    mt = mt.filter_rows(mt.locus.in_autosome())
    #remove samples which have excessive SNPs
    high_snps_file = annotdir + "high_snp_samples.txt"
    high_snp_ht = hl.import_table(high_snps_file, delimiter="\t").key_by('s')
    mt = mt.filter_cols(hl.is_defined(high_snp_ht[mt.s]), keep=False)

    #filter MT by depth/gq/vaf
    if min_dp > 0 or min_gq > 0 or min_vaf > 0:
        vaf = mt.AD[1] / hl.sum(mt.AD)
        print("Filtering input MT by depth: DP=" + str(min_dp) + ", genotype quality: GQ=" + str(min_gq) + ", VAF: VAF=" + str(min_vaf))
        filter_condition = ( (mt.GT.is_het() & (vaf > min_vaf) & (mt.DP > min_dp) & (mt.GQ > min_gq)) | 
                   (mt.GT.is_hom_ref() & (mt.DP > min_dp) & (mt.GQ > min_gq)) |
                   (mt.GT.is_hom_var() & (mt.DP > min_dp) & (mt.GQ > min_gq)))
        fraction_filtered = mt.aggregate_entries(hl.agg.fraction(~filter_condition))
        print(f'Filtering {fraction_filtered * 100:.2f}% entries out of downstream analysis.')
        mt = mt.filter_entries(filter_condition)

    mt = hl.sample_qc(mt)
    print("Writing annotated MT to file")
    mt.write(mt_sample_qc_file, overwrite=True)
    # also write just the cols as an HT
    mt.cols().write(ht_sample_qc_file,  overwrite=True)


def main():
    #set up
    inputs = parse_config()
    mtdir = inputs['matrixtables_lustre_dir']
    annotdir = inputs['annotation_lustre_dir']
    resourcesdir = inputs['resource_dir']

    #initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    #annotate mt with runid and pop
    raw_mt_file = mtdir + "gatk_unprocessed.mt"
    pop_ht_file = mtdir + "pop_assignments.ht"
    runid_file = resourcesdir + "sequencing_batches.txt"
    annotated_mt_file = mtdir + "gatk_unprocessed_with_pop_and_runid.mt"
    annotate_mt(raw_mt_file, pop_ht_file, runid_file, annotated_mt_file)

    #run sample QC without filtering entries
    mt_sample_qc_file = mtdir + "non_stratified_sample_qc.mt"
    ht_sample_qc_file = mtdir + "non_stratified_sample_qc_cols.ht"
    run_sample_qc(annotated_mt_file, mt_sample_qc_file, ht_sample_qc_file, annotdir)

    #run sample QC with sringent filters
    mt_sample_qc_file_s = mtdir + "non_stratified_sample_qc_dp20_gq20_vaf0_25.mt"
    ht_sample_qc_file_s = mtdir + "non_stratified_sample_qc_cols_dp20_gq20_vaf0_25.ht"
    run_sample_qc(annotated_mt_file, mt_sample_qc_file_s, ht_sample_qc_file_s, annotdir, 20, 20, 0.25)

    #run sample QC with more relaxed filters
    mt_sample_qc_file_r = mtdir + "non_stratified_sample_qc_dp10_gq10_vaf0_15.mt"
    ht_sample_qc_file_r = mtdir + "non_stratified_sample_qc_cols_dp10_gq10_vaf0_15.ht"
    run_sample_qc(annotated_mt_file, mt_sample_qc_file_r, ht_sample_qc_file_r, annotdir, 10, 10, 0.15)


if __name__ == '__main__':
    main()