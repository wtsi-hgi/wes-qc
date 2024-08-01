# perform hail sample QC stratified by superpopulation and identify outliers
import hail as hl
import pyspark
import os
from utils.utils import parse_config
from gnomad.sample_qc.filtering import compute_stratified_metrics_filter


def annotate_mt(raw_mt_file: str, pop_ht_file: str, annotated_mt_file: str, pop_pandas_file: str):
    '''
    Annotate mt with superpopulation and sequencing runid
    :param str raw_mt_file: raw mt file
    :param str pop_ht_file: file with population annotations
    :param str runid_file: metadata file to annotate run ids
    :param str annotated_mt_file: annotated mt file
    :param str pop_pandas_file: tsv from pandas df - EGA and pop
    '''
    mt = hl.read_matrix_table(raw_mt_file)
    #pop_ht = hl.read_table(pop_ht_file)
    pop_ht = hl.import_table(pop_pandas_file, delimiter="\t").key_by('s')
    mt = mt.annotate_cols(assigned_pop=pop_ht[mt.s].pop)
    mt.write(annotated_mt_file, overwrite=True)


def stratified_sample_qc(annotated_mt_file: str, mt_qc_outfile: str, ht_qc_cols_outfile: str, qc_filter_file: str, annotdir: str):
    '''
    Run sample QC and stratify by population
    param str annotated_mt_file: population and run id annotated MT file
    param str mt_qc_outfile: sample QC MT file
    param str ht_qccols_outfile: sample QC columns HT file
    param str qc_filter_file: output file for stratified sample QC HT
    param str annotdir: output directory for annotations
    '''
    mt = hl.read_matrix_table(annotated_mt_file)
    # filter to autosomes only
    mt = mt.filter_rows(mt.locus.in_autosome())

    # filter MT by depth/gq/vaf

    # TODO: move to config
    min_dp = 20
    min_gq = 20
    min_vaf = 0.25
    if min_dp > 0 or min_gq > 0 or min_vaf > 0:
        vaf = mt.AD[1] / hl.sum(mt.AD)
        print("Filtering input MT by depth: DP=" + str(min_dp) +
              ", genotype quality: GQ=" + str(min_gq) + ", VAF: VAF=" + str(min_vaf))
        filter_condition = ((mt.GT.is_het() & (vaf > min_vaf) & (mt.DP > min_dp) & (mt.GQ > min_gq)) |
                            (mt.GT.is_hom_ref() & (mt.DP > min_dp) & (mt.GQ > min_gq)) |
                            (mt.GT.is_hom_var() & (mt.DP > min_dp) & (mt.GQ > min_gq)))
        fraction_filtered = mt.aggregate_entries(
            hl.agg.fraction(~filter_condition))
        print(
            f'Filtering {fraction_filtered * 100:.2f}% entries out of downstream analysis.')
        mt = mt.filter_entries(filter_condition)

    # run sample QC
    mt_with_sampleqc = hl.sample_qc(mt, name='sample_qc')
    # annotate with heterozygosity rate
    mt_with_sampleqc = mt_with_sampleqc.annotate_cols(sample_qc=mt_with_sampleqc.sample_qc.annotate(
        heterozygosity_rate=mt_with_sampleqc.sample_qc.n_het/mt_with_sampleqc.sample_qc.n_called))
    mt_with_sampleqc.write(mt_qc_outfile, overwrite=True)
    mt_with_sampleqc.cols().write(ht_qc_cols_outfile,  overwrite=True)
    # stratify by pop
    pop_ht = hl.read_table(ht_qc_cols_outfile)
    qc_metrics = ['heterozygosity_rate', 'n_snp', 'r_ti_tv', 'n_transition', 'n_transversion',
                  'r_insertion_deletion', 'n_insertion', 'n_deletion', 'r_het_hom_var']

    pop_ht_afr = pop_ht.filter(pop_ht.assigned_pop == 'AFR')
    pop_ht_nonafr = pop_ht.filter(pop_ht.assigned_pop != 'AFR')

    pop_filter_ht_afr = compute_stratified_metrics_filter(
        pop_ht_afr,
        qc_metrics={
            metric: pop_ht_afr.sample_qc[metric] for metric in qc_metrics
        },
        metric_threshold={'heterozygosity_rate': (100, 4),
                          'r_het_hom_var': (100, 5)},

        strata={"qc_pop": pop_ht_afr.assigned_pop},
    )

    pop_filter_ht_nonafr = compute_stratified_metrics_filter(
        pop_ht_nonafr,
        qc_metrics={
            metric: pop_ht_nonafr.sample_qc[metric] for metric in qc_metrics
        },
        metric_threshold={'heterozygosity_rate': (100, 4),
                          'r_het_hom_var': (100, 4)},

        strata={"qc_pop": pop_ht_nonafr.assigned_pop},
    )

    pop_filter_ht = pop_filter_ht_afr.union(pop_filter_ht_nonafr)
    globals1 = hl.eval(pop_filter_ht_nonafr.globals.qc_metrics_stats)
    globals2 = hl.eval(pop_filter_ht_afr.globals.qc_metrics_stats)
    globals1.update(globals2)

    pop_ht = pop_ht.annotate_globals(qc_metrics_stats=globals1)
    pop_ht = pop_ht.annotate(**pop_filter_ht[pop_ht.key]).persist()
    checkpoint = pop_ht.aggregate(hl.agg.count_where(
        hl.len(pop_ht.qc_metrics_filters) == 0))
    print(f'{checkpoint} exome samples found passing pop filtering')
    pop_ht.write(qc_filter_file, overwrite=True)

    output_text_file = os.path.join(annotdir, "sample_qc_by_pop.tsv.bgz")
    pop_ht.export(output_text_file, delimiter="\t")

    output_globals_json = os.path.join(annotdir, "sample_qc_by_pop.globals.json")
    pop_ht.globals.export(output_globals_json)


def main():
    # set up
    inputs = parse_config()
    mtdir = inputs['matrixtables_lustre_dir']
    annotdir = inputs['annotation_lustre_dir']
    resourcesdir = inputs['resource_dir']
    #mtdir2 = inputs['load_matrixtables_lustre_dir']

    # initialise hail
    tmp_dir = inputs['tmp_dir']
    sc = pyspark.SparkContext.getOrCreate()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38", idempotent=True)

    # annotate mt with runid and pop
    raw_mt_file = os.path.join(mtdir, "gatk_unprocessed.mt")
    pop_ht_file = os.path.join(mtdir, "pop_assignments.ht")
    annotated_mt_file = os.path.join(mtdir, "gatk_unprocessed_with_pop.mt")
    pop_pandas_file = os.path.join(mtdir, "pop_assignments.tsv")
    annotate_mt(raw_mt_file, pop_ht_file, annotated_mt_file, pop_pandas_file)

    # run sample QC and stratify by population
    mt_qc_outfile = os.path.join(mtdir, "mt_pops_sampleqc.mt")
    ht_qc_cols_outfile = os.path.join(mtdir, "mt_pops_sampleqc.ht")
    qc_filter_file = os.path.join(mtdir, "mt_pops_QC_filters.ht")
    stratified_sample_qc(annotated_mt_file, mt_qc_outfile,
                         ht_qc_cols_outfile, qc_filter_file, annotdir)


if __name__ == '__main__':
    main()
