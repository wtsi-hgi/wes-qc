# perform hail sample QC on raw matrixtable and annotate with superpops and fail count from stratified qc
import resource
import hail as hl
import pyspark
from wes_qc.utils.utils import parse_config

def run_sample_qc(raw_mt_file, mt_sample_qc_file, ht_sample_qc_file, metadatafile, pop_ht_file, annotdir):
    '''
    Run sample QC on raw variants
    :param str raw_mt_file: raw mt file 
    :param str mt_sample_qc_file: output MT
    :param str ht_sample_qc_file: output HT, cols only
    :param str metadatafile: metadata file to annotate batches
    :param str pop_ht_file: file with popuklation annotations
    :param str annotdir: annotations dir
    '''
    mt = hl.read_matrix_table(raw_mt_file)
    #filter to autosomes only
    mt = mt.filter_rows(mt.locus.in_autosome())
    #remove samples which have excessive SNPs
    high_snps_file = annotdir + "high_snp_samples.txt"
    high_snp_ht = hl.import_table(high_snps_file, delimiter="\t").key_by('s')
    mt = mt.filter_cols(hl.is_defined(high_snp_ht[mt.s]), keep=False)
    #annotate with opoulation, sequencing batch and location 
    metadata_ht = hl.import_table(metadatafile, delimiter="\t").key_by('ega')
    mt = mt.annotate_cols(batch=metadata_ht[mt.s]['batch.y'])
    seq_expr = (hl.case()
                   .when(mt.s.startswith('EGAN'), 'Sanger')
                   .when(mt.s.startswith('Z'), 'Bristol')
                   .default("")
                  )
    mt = mt.annotate_cols(sequencing_location = seq_expr).key_cols_by('s')
    pop_ht = hl.read_table(pop_ht_file)
    mt = mt.annotate_cols(assigned_pop=pop_ht[mt.s].pop)

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

    #run sample QC
    raw_mt_file = mtdir + "gatk_unprocessed.mt"
    mt_sample_qc_file = mtdir + "non_stratified_sample_qc.mt"
    ht_sample_qc_file = mtdir + "non_stratified_sample_qc_cols.ht"
    metadatafile = resourcesdir + "all_samples_with_proceed_and_seq_info_and_warehouse_info_egas_from_mlwh.txt"
    pop_ht_file = mtdir + "pop_assignments.ht"
    run_sample_qc(raw_mt_file, mt_sample_qc_file, ht_sample_qc_file, metadatafile, pop_ht_file, annotdir)


if __name__ == '__main__':
    main()