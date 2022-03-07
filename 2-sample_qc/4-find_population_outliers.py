# perform hail sample QC stratified by superpopulation and identify outliers
import hail as hl
import pyspark
from wes_qc.utils.utils import parse_config
from gnomad.sample_qc.filtering import compute_stratified_metrics_filter


def annotate_pops(pruned_mt_file: str, pop_ht_file: str, pops_annot_mt_file: str):
    '''
    Annotate LD pruned matrixtable with predicted populations
    :param str pruned_mt_file: ld pruned MT file
    :param str pop_ht_file: predicted population HT file
    :param str pops_annot_mt_file: population annotated MT file
    '''
    mt = hl.read_matrix_table(pruned_mt_file)
    pop_ht = hl.read_table(pop_ht_file)
    mt = mt.annotate_cols(assigned_pop=pop_ht[mt.s].pop)
    mt.write(pops_annot_mt_file, overwrite=True)


def stratified_sample_qc(pops_annot_mt_file: str, mt_qc_outfile: str, ht_qc_cols_outfile: str, qc_filter_file: str, annotdir: str):
    '''
    Run sample QC and stratify by population
    param str pops_annot_mt_file: population annotated MT file
    param str mt_qc_outfile: sample QC MT file
    param str ht_qccols_outfile: sample QC columns HT file
    param str qc_filter_file: output file for stratified sample QC HT
    aram str annotdir: output directory for annotations
    '''
    mt = hl.read_matrix_table(pops_annot_mt_file)
    #run sample QC
    mt_with_sampleqc = hl.sample_qc(mt, name='sample_qc')
    #annotate with heterozygosity rate
    mt_with_sampleqc = mt_with_sampleqc.annotate_cols(sample_qc=mt_with_sampleqc.sample_qc.annotate(
        heterozygosity_rate=mt_with_sampleqc.sample_qc.n_het/mt_with_sampleqc.sample_qc.n_called))
    mt_with_sampleqc.write(mt_qc_outfile, overwrite=True)
    mt_with_sampleqc.cols().write(ht_qc_cols_outfile,  overwrite=True)
    #stratify by pop
    pop_ht = hl.read_table(ht_qc_cols_outfile)
    qc_metrics = ['heterozygosity_rate', 'n_snp', 'r_ti_tv',
                  'r_insertion_deletion', 'n_insertion', 'n_deletion', 'r_het_hom_var']

    pop_filter_ht = compute_stratified_metrics_filter(
        pop_ht, 
        qc_metrics={
            metric: pop_ht.sample_qc[metric] for metric in qc_metrics
        },
        strata={"qc_pop": pop_ht.assigned_pop},
        )
    
    pop_ht = pop_ht.annotate_globals(**pop_filter_ht.index_globals())
    pop_ht = pop_ht.annotate(**pop_filter_ht[pop_ht.key]).persist()
    checkpoint = pop_ht.aggregate(hl.agg.count_where(
        hl.len(pop_ht.qc_metrics_filters) == 0))
    print(f'{checkpoint} exome samples found passing pop filtering')
    pop_ht.write(qc_filter_file, overwrite=True)

    output_text_file = annotdir + "sample_qc_by_pop.tsv.bgz"
    pop_ht.export(output_text_file, delimiter="\t")


def main():
    #set up
    inputs = parse_config()
    mtdir = inputs['matrixtables_lustre_dir']
    annotdir = inputs['annotation_lustre_dir']

    #initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    #add population labels to ld-pruned biallelic autosomal SNPs
    pruned_mt_file = mtdir + "mt_ldpruned.mt"
    pop_ht_file = mtdir + "pop_assignments.ht"
    pops_annot_mt_file = mtdir + "mt_ldpruned_pops_annotated.mt"
    annotate_pops(pruned_mt_file, pop_ht_file, pops_annot_mt_file)

    #run sample QC and stratify by population
    mt_qc_outfile = mtdir + "mt_pops_sampleqc.mt"
    ht_qc_cols_outfile = mtdir + "mt_pops_sampleqc.ht"
    qc_filter_file = mtdir + "mt_pops_QC_filters.ht"
    stratified_sample_qc(pops_annot_mt_file, mt_qc_outfile, ht_qc_cols_outfile, qc_filter_file, annotdir)


if __name__ == '__main__':
    main()
