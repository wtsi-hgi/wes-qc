# median counts of variants of each consequence per sample plus transmitted/untransmitted 
# ratio for synonymous singletons
import hail as hl
import pyspark
import argparse
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
    synon_mt = mt_in.filter_rows(mt_in.info.consequence == 'synonymous_variant')
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
    #split by consequence
    median_count_for_cq(snv_mt, ['synonymous_variant', 'misense_variant'])
    # synonymous_mt = snv_mt.filter_rows(snv_mt.info.consequence == 'synonymous_variant')
    # missense_mt = snv_mt.filter_rows(snv_mt.info.consequence == 'missense_variant')
    # nonsense_mt = snv_mt.filter_rows(snv_mt.info.consequence == 'stop_gained')
    # splice_acc_donor_mt = snv_mt.filter_rows((snv_mt.info.consequence == 'splice_acceptor_variant') | (snv_mt.info.consequence == 'splice_donor_variant') )
    # frameshift_mt = indel_mt.filter_rows(indel_mt.info.consequence == 'frameshift_variant')
    # inframe_mt = indel_mt.filter_rows((indel_mt.info.consequence == 'inframe_deletion') | (indel_mt.info.consequence == 'inframe_insertion') )

    # coding_snvs = ['synonymous_variant', 'missense_variant', 'stop_gained','splice_acceptor_variant', 'splice_donor_variant', 'sart_lost', 'stop_lost']
    # coding_indels = ['frameshift_variant', 'inframe_deletion', 'inframe_insertion']
    # coding_snv_mt = snv_mt.filter_rows(snv_mt.info.consequence in coding_snvs)
    # coding_indel_mt = indel_mt.filter_rows(indel_mt.info.consequence in coding_indels)

def median_count_for_cq(mt_in: hl.MatrixTable, cqs: list):
    '''
    Get median counts per sample for a list of consequences
    :param hl.MatrixTable mt_in: Input MatrixTable
    :param list cqs: List fof consequences
    '''
   # mt = mt_in.filter_rows(mt_in.info.consequence in cqs)
    mt = mt_in.filter_rows(hl.literal(mt_in.info.consequence).contains(cqs))
    x = mt.aggregate_rows(hl.agg.counter(mt.info.consequence))
    x = dict(x)
    print(x)
    exit(0)
    mt_rare = mt.filter_rows(mt.gnomad_AC < 5)
    mt = hl.sample_qc(mt)
    mt_rare = hl.sample_qc(mt_rare)
    sampleqc_ht = mt.cols()
    sampleqc_rare_ht = mt_rare.cols()
    print(hl.median(sampleqc_ht.sample_qc.n_non_ref))
    print(hl.median(sampleqc_rare_ht.sample_qc.n_non_ref))




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