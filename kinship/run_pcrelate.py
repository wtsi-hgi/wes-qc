#apply RF model for variant QC
import hail as hl
import pyspark
import argparse
from wes_qc.utils.utils import parse_config


def filter_mt(mtfile: str, mtfile_for_pcrelate: str):
    '''
    :param str mtfile: Input matrixtable file
    :param str mtfile_for_pcrelate: Output matrixtable file
    Filter mtfile before running PC relate - split multi, apply hard filters used for variant QC, autosomes
    '''
    mt = hl.read_matrix_table(mtfile)
    # filter to autosomes only
    mt = mt.filter_rows(mt.locus.in_autosome())
    #split mutli
    mt = split_multi_hts(mt)
    # filter MT 
    filter_condition = ((mt.info.QD < 2) | (mt.info.FS > 60) | (mt.info.MQ < 30))
    mt = mt.filter_rows(filter_condition, keep=False)

    mt.write(mtfile_for_pcrelate, overwrite=True)


def run_pcrelate(mtfile_for_pcrelate, pcrelate_htoutfile, pcrelate_textfile):
    '''
    Run PCrelate and save output table
    :param str mtfile_for_pcrelate: input mtfile
    :param str pcrelate_hroutfile: output hail table file
    :paran str pcrelate_textfile: output textfile
    '''
    mt = hl.read_matrix_table(mtfile_for_pcrelate)
    pcrht = hl.pc_relate(mt.GT, 0.01, k=10)
    pcrht.write(pcrelate_htoutfile, overwrite=True)
    #text file
    pcrht.export(pcrelate_textfile, delimiter="\t")


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

    mtfile = mtdir + "gatk_unprocessed_with_pop_and_runid.mt"#use an unfiltered mt, pre-sample QC
    mtfile_for_pcrelate = mtdir + "mt_for_pcrelate.mt"
    pcrelate_outfile = mtdir + "pcrelate.ht"
    pcrelate_textfile = annotdir + "pcrelate.tsv.bgz"

    filter_mt(mtfile, mtfile_for_pcrelate)
    run_pcrelate(mtfile_for_pcrelate, pcrelate_outfile, pcrelate_textfile)


if __name__ == '__main__':
    main() 