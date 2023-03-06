#apply hard filters to genotypes
import hail as hl
import pyspark
import argparse
from utils.utils import parse_config

def get_options():
    '''
    Get options from the command line
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("--dp", type = int, help="DP cut off")
    parser.add_argument("--gq", type = int, help="GQ cut off")
    parser.add_argument("--ab", type = float, help="Heterozygous allele balance cut off")
    args = parser.parse_args()
    if not args.dp and args.gq and args.ab:
        print("--dp, --gq and --ab must be specified")
        exit(1)

    return args


def filter_mt(mtfile: str, dp: int, gq: int, ab: float, mtfile_filtered: str):
    '''
    Filter with hard genotype filters
    :param str mtfile: Input mtfile
    :param int dp: Minimum DP
    :param int gq: Minimum GQ
    :param float ab: Minimum allele balance for hets
    :param str mtfile_filtered: Output mtfile
    '''
    mt = hl.read_matrix_table(mtfile)
    prev_count = mt.entries().count()
    var_count = mt.count_rows()

    filter_condition = (
        (mt.GT.is_het() & (mt.HetAB < ab)) | 
        (mt.DP < dp) |
        (mt.GQ < gq)
    )
    mt = mt.annotate_entries(
        hard_filters = hl.if_else(filter_condition, 'Fail', 'Pass')
    )
    mt = mt.filter_entries(mt.hard_filters == 'Pass')

    
    # mt = mt.filter_entries(
    #     (mt.GT.is_het() & (mt.HetAB >= ab)) | 
    #     (mt.DP > dp) | 
    #     (mt.GQ > gq)
    # )
    #remove unused rows
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(mt.variant_qc.n_non_ref == 0, keep = False)
    new_count = mt.entries().count()
    new_var_count = mt.count_rows()

    mt.write(mtfile_filtered, overwrite = True)
    print(str(prev_count) + " entries prior to filtering, " + str(new_count) + " entries after filtering")
    print(str(var_count) + " variants prior to filtering, " + str(new_var_count) + " variants after filtering")


def main():
    #set up
    args = get_options()
    inputs = parse_config()
    mtdir = inputs['matrixtables_lustre_dir']

    # initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    mtfile = mtdir + "mt_after_var_qc.mt"
    mtfile_filtered = mtdir + "mt_after_var_qc_hard_filter_gt.mt"
    filter_mt(mtfile, args.dp, args.gq, args.ab, mtfile_filtered)


if __name__ == '__main__':
    main()