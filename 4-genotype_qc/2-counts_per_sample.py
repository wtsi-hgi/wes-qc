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
    synon_mt = mt_in.filter_rows(synon_mt.info.consequence == 'synonymous_variant')
    tdt_ht = hl.transmission_disequilibrium_test(mt_in, pedigree)

    trans_sing = tdt_ht.filter((tdt_ht.t == 1) & (tdt_ht.u == 0))
    trans = trans_sing.count()

    untrans_sing = tdt_ht.filter((tdt_ht.t == 0) & (tdt_ht.u == 1))
    untrans = untrans_sing.count()
    ratio = trans/untrans

    print("There are " + str(trans) + " transmitted synonymous singletons and " + str(untrans) + " untransmitted synonymous singletons")
    print("The ratio of transmitted/unstransmitted singletons is " +'{0:.4f}'.format(ratio) )


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



if __name__ == '__main__':
    main()