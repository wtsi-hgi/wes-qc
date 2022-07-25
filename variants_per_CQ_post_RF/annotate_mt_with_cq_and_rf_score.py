#
import hail as hl
import pyspark
import argparse
from wes_qc.utils.utils import parse_config


def get_options():
    '''
    Get options from the command line
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("--runhash", help="RF training run hash")
    args = parser.parse_args()
    if not args.runhash:
        print("--runhash must be specified")
        exit(1)

    return args


def annotate_mt_with_consequence(mtfile: str, cqfile: str, cq_annotated_mtfile: str):
    '''
    Annotate matrixtable with consequences
    :param str mtfile: input matrixtable file
    :param str cqfile: input VEP consequence file
    :param str cq_annotated_mtfile: output consequence annotated mtfile
    '''
    mt = hl.read_matrix_table(mtfile)

    cq_ht=hl.import_table(cqfile,types={'f0':'str','f1':'int32', 'f2':'str','f3':'str','f4':'str', 'f5':'str'}, no_header=True)
    cq_ht=cq_ht.annotate(chr=cq_ht.f0)
    cq_ht=cq_ht.annotate(pos=cq_ht.f1)
    cq_ht=cq_ht.annotate(rs=cq_ht.f2)
    cq_ht=cq_ht.annotate(ref=cq_ht.f3)
    cq_ht=cq_ht.annotate(alt=cq_ht.f4)
    cq_ht=cq_ht.annotate(consequence=cq_ht.f5)
    cq_ht = cq_ht.key_by(
    locus=hl.locus(cq_ht.chr, cq_ht.pos), alleles=[cq_ht.ref,cq_ht.alt])
    cq_ht=cq_ht.drop(cq_ht.f0,cq_ht.f1,cq_ht.f2,cq_ht.f3,cq_ht.f4,cq_ht.chr,cq_ht.pos,cq_ht.ref,cq_ht.alt)
    cq_ht = cq_ht.key_by(cq_ht.locus, cq_ht.alleles)

    mt = mt.annotate_rows(
        info=mt.info.annotate(
            consequence=cq_ht[mt.row_key].consequence)
    )
    mt.write(cq_annotated_mtfile, overwrite = True)


def annotate_mt_with_rf_score(cq_annotated_mtfile: str, rf_htfile: str, rf_annotated_mtfile: str):
    '''
    Annotate matrixtable with RF score and bin
    :param str cq_annotated_mtfile: consequence annotated mtfile
    :param str rf_htfile: random forest hail table file
    :param str rf_annotated_mtfile: random forest score annotated mtfile
    '''
    mt = hl.read_matrix_table(cq_annotated_mtfile)
    rf_ht = hl.read_table(rf_htfile)

    mt = mt.annotate_rows(
        info=mt.info.annotate(
            rf_probability=rf_ht[mt.row_key].rf_probability['TP'])
    )
    
    mt.write(rf_annotated_mtfile, overwrite=True)


def main():
    # set up
    args = get_options()
    inputs = parse_config()
    rf_dir = inputs['var_qc_rf_dir']
    mtdir = inputs['matrixtables_lustre_dir']
    resourcedir = inputs['resource_dir']

    # initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    mtfile = mtdir + "mt_varqc_splitmulti.mt"#matrix table with multiallelics split
    cqfile = resourcedir + "all_consequences.txt"
    cq_annotated_mtfile = mtdir + "mt_varqc_splitmulti_with_cq.mt"
    annotate_mt_with_consequence(mtfile, cqfile, cq_annotated_mtfile)

    rf_annotated_mtfile = mtdir + "mt_varqc_splitmulti_with_cq_and_rf_scores_" + args.run_hash + ".mt"
    rf_htfile = rf_dir + args.runhash + "/rf_result_final_for_ranking.ht"
    annotate_mt_with_rf_score(cq_annotated_mtfile, rf_htfile, rf_annotated_mtfile)


if __name__ == '__main__':
    main()