#find population outliers based on sequencing location and based on sequencing location + superpopulation
#uses sample QC hail table from 4-find_population_outliers.py
import hail as hl
import pyspark
from wes_qc.utils.utils import parse_config
from gnomad.sample_qc.filtering import compute_stratified_metrics_filter


def stratify_sample_qc(ht_qc_cols_outfile, qc_filter_file_seq_loc, qc_filter_file_seq_loc_pop, qc_filter_seq_text_file, qc_filter_seq_pop_text_file):
    '''
    Stratify sample QC by sequencing location, and by sequencing location + pop
    :param str ht_qc_cols_outfile: hail table file with sample QC
    :param str qc_filter_file_seq_loc: outfile for sample QC stratified by sequencing location
    :param str qc_filter_file_seq_loc_pop: outfile for sample QC stratified by sequencing location and superpopulation
    :param str qc_filter_seq_text_file: text outfile for sample QC stratified by sequencing location
    :param str qc_filter_seq_pop_text_file: text outfile for sample QC stratified by sequencing location and superpopulation
    '''
    qc_ht = hl.read_table(ht_qc_cols_outfile)
    #annotate sample qc ht with strata based on sequencing location and superpop
    strat_expr = (hl.case()
                   .when(qc_ht.sequencing_location == 'Bristol', 'Broad')
                   .when(qc_ht.sequencing_location == 'Sanger', qc_ht.assigned_pop)
                   .default("")
                  )
    qc_ht = qc_ht.annotate(strat = strat_expr).key_by('s')

    qc_metrics = ['heterozygosity_rate', 'n_snp', 'r_ti_tv', 'n_transition', 'n_transversion',
                  'r_insertion_deletion', 'n_insertion', 'n_deletion', 'r_het_hom_var']

    #stratify based on sequencing location only
    pop_filter_ht_loc = compute_stratified_metrics_filter(
        qc_ht,
        qc_metrics={
            metric: qc_ht.sample_qc[metric] for metric in qc_metrics
        },
        strata={"qc_seq_loc": qc_ht.sequencing_location},
    )

    qc_ht2 = qc_ht
    qc_ht2 = qc_ht2.annotate_globals(**pop_filter_ht_loc.index_globals())
    qc_ht2 = qc_ht2.annotate(**pop_filter_ht_loc[qc_ht2.key]).persist()
    checkpoint = qc_ht2.aggregate(hl.agg.count_where(
        hl.len(qc_ht2.qc_metrics_filters) == 0))
    print(f'{checkpoint} exome samples found passing pop filtering stratified by sequencing location only')
    qc_ht2.write(qc_filter_file_seq_loc, overwrite=True)
    qc_ht2.export(qc_filter_seq_text_file, delimiter="\t")

    #stratify based on sequencing location and superpopulation
    pop_filter_ht_loc_pop = compute_stratified_metrics_filter(
        qc_ht,
        qc_metrics={
            metric: qc_ht.sample_qc[metric] for metric in qc_metrics
        },
        strata={"qc_seq_loc_pop": qc_ht.strat},
    )
    qc_ht3 = qc_ht
    qc_ht3 = qc_ht3.annotate_globals(**pop_filter_ht_loc_pop.index_globals())
    qc_ht3 = qc_ht3.annotate(**pop_filter_ht_loc_pop[qc_ht3.key]).persist()
    checkpoint = qc_ht3.aggregate(hl.agg.count_where(
        hl.len(qc_ht3.qc_metrics_filters) == 0))
    print(f'{checkpoint} exome samples found passing pop filtering stratified by sequencing location and superpopulation')
    qc_ht3.write(qc_filter_file_seq_loc_pop, overwrite=True)
    qc_ht3.export(qc_filter_seq_pop_text_file, delimiter="\t")
    

def main():
    # set up
    inputs = parse_config()
    mtdir = inputs['matrixtables_lustre_dir']
    annotdir = inputs['annotation_lustre_dir']
    resourcesdir = inputs['resource_dir']

    # initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    #stratify sample QC
    ht_qc_cols_outfile = mtdir + "mt_pops_sampleqc.ht"
    qc_filter_file_seq_loc = mtdir + "mt_pops_QC_filters_sequencing_location.ht"
    qc_filter_file_seq_loc_pop = mtdir + "mt_pops_QC_filters_sequencing_location_and_superpop.ht"
    qc_filter_seq_text_file = annotdir + "sample_qc_by_seq_location.tsv.bgz"
    qc_filter_seq_pop_text_file = annotdir + "sample_qc_by_seq_location_and_superpop.tsv.bgz"
    stratify_sample_qc(ht_qc_cols_outfile, qc_filter_file_seq_loc, qc_filter_file_seq_loc_pop, qc_filter_seq_text_file, qc_filter_seq_pop_text_file)


if __name__ == '__main__':
    main()