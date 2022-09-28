# median counts of variants of each consequence per sample plus transmitted/untransmitted 
# ratio for synonymous singletons
import hail as hl
import pyspark
import argparse
import pandas as pd
from wes_qc.utils.utils import parse_config


def annotate_gnomad(mt_in: hl.MatrixTable, gnomad_htfile: str) -> hl.MatrixTable:
    '''
    Annotate matrixtable with AC and AF from gnomad
    :param hl.MatrixTable mt_in: Input MatrixTable
    :param str gnomad_htfile: gnomAD Hail table file
    :return: Annotated MatrixTable
    '''
    gnomad_ht = hl.read_table(gnomad_htfile)
    #ac = gnomadht.freq[0].AC
    #af = gnomadht.freq[0].AF
    mt_in = mt_in.annotate_rows(gnomad_AF=gnomad_ht[mt_in.row_key].freq[0].AF)
    mt_in = mt_in.annotate_rows(gnomad_AC=gnomad_ht[mt_in.row_key].freq[0].AC)

    return mt_in
    

def get_trans_untrans_synon_singleton_counts(mt_in: hl.MatrixTable, pedfile: str):
    '''
    Get transmitted/untransmitted counts and ratio for synonymous singletons
    :param hl.MatrixTable mt_in: Input MatrixTable
    :param str pedfile: Path to pedfile
    '''
    pedigree = hl.Pedigree.read(pedfile)
    #get trios only then filter to trioAC == 1 or 2
    trio_sample_ht = hl.import_fam(pedfile)
    sample_list = trio_sample_ht.id.collect() + trio_sample_ht.pat_id.collect() + trio_sample_ht.mat_id.collect()
    mt2 = mt_in.filter_cols(hl.set(sample_list).contains(mt_in.s))
    mt2 = hl.variant_qc(mt2, name = 'varqc_trios')
    mt2 = mt2.filter_rows(mt2.varqc_trios.AC[1] <= 2)

    synon_mt = mt2.filter_rows(mt2.info.consequence == 'synonymous_variant')
    tdt_ht = hl.transmission_disequilibrium_test(synon_mt, pedigree)

    trans_sing = tdt_ht.filter((tdt_ht.t == 1) & (tdt_ht.u == 0))
    trans = trans_sing.count()

    untrans_sing = tdt_ht.filter((tdt_ht.t == 0) & (tdt_ht.u == 1))
    untrans = untrans_sing.count()
    ratio = trans/untrans

    print("There are " + str(trans) + " transmitted synonymous singletons and " + str(untrans) + " untransmitted synonymous singletons")
    print("The ratio of transmitted/unstransmitted singletons is " +'{0:.4f}'.format(ratio) )


def get_counts_per_cq(mt_in: hl.MatrixTable):
    '''
    Get median counts of each consequence per sample
    :param hl.MatrixTable mt_in: Input MatrixTable
    '''
    #split mt by snvs and indels
    snv_mt = mt_in.filter_rows(hl.is_snp(mt_in.alleles[0], mt_in.alleles[1]))
    indel_mt = mt_in.filter_rows(hl.is_indel(mt_in.alleles[0], mt_in.alleles[1]))

    #get median numbers of variants by consequence
    synonymous_counts =  median_count_for_cq(snv_mt, ['synonymous_variant'])
    missense_counts =  median_count_for_cq(snv_mt, ['missense_variant'])
    nonsense_counts =  median_count_for_cq(snv_mt, ['stop_gained'])
    splice_acc_don_counts = median_count_for_cq(snv_mt, ['splice_acceptor_variant', 'splice_donor_variant'])

    frameshift_counts =  median_count_for_cq(indel_mt, ['frameshift_variant'])
    inframe_indel_counts = median_count_for_cq(indel_mt, ['inframe_deletion', 'inframe_insertion'])

    coding_snv_counts = median_count_for_cq(snv_mt, ['synonymous_variant', 'missense_variant', 'stop_gained','splice_acceptor_variant', 'splice_donor_variant', 'sart_lost', 'stop_lost'])
    coding_indel_counts = median_count_for_cq(indel_mt, ['frameshift_variant', 'inframe_deletion', 'inframe_insertion'])

    print("Synonymous: Total " + str(synonymous_counts[0]) + " rare " + str(synonymous_counts[1]))
    print("Missense: Total " + str(missense_counts[0]) + " rare " + str(missense_counts[1]))
    print("Nonsense: Total " + str(nonsense_counts[0]) + " rare " + str(nonsense_counts[1]))
    print("Splicing: Total " + str(splice_acc_don_counts[0]) + " rare " + str(splice_acc_don_counts[1]))
    print("Frameshift: Total " + str(frameshift_counts[0]) + " rare " + str(frameshift_counts[1]))
    print("In-frame indel: Total " + str(inframe_indel_counts[0]) + " rare " + str(inframe_indel_counts[1]))
    print("Coding SNV: Total " + str(coding_snv_counts[0]) + " rare " + str(coding_snv_counts[1]))
    print("Coding indel: Total " + str(coding_indel_counts[0]) + " rare " + str(coding_indel_counts[1]))


def median_count_for_cq(mt_in: hl.MatrixTable, cqs: list) -> tuple:
    '''
    Get median counts per sample for a list of consequences
    :param hl.MatrixTable mt_in: Input MatrixTable
    :param list cqs: List fof consequences
    :return: tuple
    '''

    mt = mt_in.filter_rows(hl.literal(cqs).contains(mt_in.info.consequence))
    mt_rare = mt.filter_rows(mt.gnomad_AC < 5)

    mt = hl.sample_qc(mt)
    mt_rare = hl.sample_qc(mt_rare)
    sampleqc_ht = mt.cols()
    sampleqc_rare_ht = mt_rare.cols()
    total_median = hl.median(sampleqc_ht.sample_qc.n_non_ref.collect()).collect()[0]
    rare_median = hl.median(sampleqc_rare_ht.sample_qc.n_non_ref.collect()).collect()[0]
    
    return total_median, rare_median

def main():
    # set up
    inputs = parse_config()
    mtdir = inputs['matrixtables_lustre_dir']
    resourcedir = inputs['resource_dir']

    # initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    mtfile = mtdir + "mt_after_var_qc_hard_filter_gt.mt"
    gnomad_htfile = resourcedir + "gnomad.exomes.r2.1.1.sites.liftover_grch38.ht"
    mt = hl.read_matrix_table(mtfile)

    mt = annotate_gnomad(mt, gnomad_htfile)

    pedfile = resourcedir + "trios.ped"
    get_trans_untrans_synon_singleton_counts(mt, pedfile)

    get_counts_per_cq(mt)



if __name__ == '__main__':
    main()