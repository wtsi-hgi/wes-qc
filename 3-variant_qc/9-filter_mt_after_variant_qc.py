# filter split matrixtable after variant QC
import hail as hl
import pyspark
import argparse
from wes_qc.utils.utils import parse_config


def get_options():
    '''
    Get options from the command line
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("--runhash", help="RF run hash")
    parser.add_argument("--snv", type = int, help="SNV threshold")
    parser.add_argument("--indel", type = int, help="Indel threshold")
    args = parser.parse_args()
    if not args.runhash and args.snv and args.indel:
        print("--runhash, --snv and --indel must be specified")
        exit(1)

    return args


def annotate_mt_with_cq_rf_score_and_bin(mtfile: str, rf_htfile: str, snv_threshold: int, indel_threshold: int, cqfile: str, filtered_mtfile: str):
    '''
    Annotate matrixtable with RF score and bin then filter SNVs and indels according to threshold
    :param str mtfile: Inputmtfile
    :param str rf_htfile: random forest hail table file
    :param int snv_threshold: bin threshold for SNVs
    :param int indel_threshold: bin threshold for indels
    :param str cqfile: File containing consequence annotation
    :param str filtered_mtfile: random forest score annotated mtfile
    '''
    mt = hl.read_matrix_table(mtfile)
    rf_ht = hl.read_table(rf_htfile)

    # keep only vars with rank_id = rank
    rf_ht = rf_ht.filter(rf_ht.rank_id == 'rank')

    # annotate mt with score and bin
    mt = mt.annotate_rows(
        info=mt.info.annotate(
            rf_score=rf_ht[mt.row_key].score)
    )
    mt = mt.annotate_rows(
        info=mt.info.annotate(
            rf_bin=rf_ht[mt.row_key].bin)
    )

    # annotate with VEP consequence
    cq_ht = hl.import_table(cqfile, types={'f0': 'str', 'f1': 'int32', 'f2': 'str',
                            'f3': 'str', 'f4': 'str', 'f5': 'str'}, no_header=True)
    cq_ht = cq_ht.annotate(chr=cq_ht.f0)
    cq_ht = cq_ht.annotate(pos=cq_ht.f1)
    cq_ht = cq_ht.annotate(rs=cq_ht.f2)
    cq_ht = cq_ht.annotate(ref=cq_ht.f3)
    cq_ht = cq_ht.annotate(alt=cq_ht.f4)
    cq_ht = cq_ht.annotate(consequence=cq_ht.f5)
    cq_ht = cq_ht.key_by(
        locus=hl.locus(cq_ht.chr, cq_ht.pos), alleles=[cq_ht.ref, cq_ht.alt])
    cq_ht = cq_ht.drop(cq_ht.f0, cq_ht.f1, cq_ht.f2, cq_ht.f3, cq_ht.f4, cq_ht.chr, cq_ht.pos, cq_ht.ref, cq_ht.alt)
    cq_ht = cq_ht.key_by(cq_ht.locus, cq_ht.alleles)
    mt = mt.annotate_rows(
        info=mt.info.annotate(
            consequence=cq_ht[mt.row_key].consequence)
    )

    # filter by SNV and indel thresholds
    mt_filtered = mt.filter_rows(((hl.is_snp(mt.alleles[0], mt.alleles[1])) & (mt.info.rf_bin <= snv_threshold)) | (
        (hl.is_indel(mt.alleles[0], mt.alleles[1])) & (mt.info.rf_bin <= indel_threshold)))

    mt_filtered.write(filtered_mtfile, overwrite = True)

    nvars = mt.count_rows()
    nvar_filtered = mt_filtered.count_rows()
    print(str(nvars) + " variants befofre filtering, " + str(nvar_filtered) + " variants after filtering")


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

    htfile = rf_dir + args.runhash + "/_gnomad_score_binning_tmp.ht"
    mtfile = mtdir + "mt_varqc_splitmulti.mt"
    mtfile_after_varqc = mtdir + "mt_after_var_qc_snv_40_indel_60.mt"
    cqfile = resourcedir + "all_consequences.txt"

    annotate_mt_with_cq_rf_score_and_bin(mtfile, htfile, args.snv, args.indel, cqfile, mtfile_after_varqc)


if __name__ == '__main__':
    main()
