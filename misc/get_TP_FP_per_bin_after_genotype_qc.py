#get TP and FP percentage after variant and genotype qc
import hail as hl
import pyspark
import argparse
from wes_qc.utils.utils import parse_config


def annotate_with_rf(mt: hl.MatrixTable, rf_htfile: str) -> hl.MatrixTable:
    '''
    Annotate MatrixTable with TP, FP, rf_bin and rf_score
    :param hl.MatrixTable mt: Input MatrixTable
    :param str rf_htfile: Random forest ht file
    :return: hl.MatriTable
    '''
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

    mt = mt.annotate_rows(TP=rf_ht[mt.row_key].tp)
    mt = mt.annotate_rows(FP=rf_ht[mt.row_key].fp)


def count_tp_fp(mt) -> tuple:
    '''
    Count total TPs and FPs in a MatrixTable
    :param hl.MatrixTable mt: Input mt
    :return: tuple of TP and FP counts
    '''
    mt_tmp = mt.filter_rows(mt.TP == True)
    tp_count = mt_tmp.count_rows()

    mt_tmp = mt.filter_rows(mt.FP == True)
    fp_count = mt_tmp.count_rows()

    return tp_count, fp_count


def get_threshold_range(threshold: int) -> list:
    '''
    Given a single threshold value convert to a range - for each threshold also look at 3 bins on either 
    side and +/-5 (unless <1 or >101)
    :param int threshold: single threshold value
    '''
    t_list = list(range(threshold-3, threshold+4))
    min_t = t_list[0]-5
    max_t = t_list[6]+5
    t_list.insert(0,min_t)
    t_list.append(max_t)
    while t_list[len(t_list)-1] > 101:
        t_list.pop()
    while(t_list[0] < 1):
        t_list.pop[0]

    return t_list


def filter_mt_count_tp_fp(mt: hl.MatrixTable, bins: list) -> dict:
    '''
    Filter mt by each rf bin in a list, then genotype hard filters and count remaining TP and P variants
    :param hl.MatrixTable mt: Input MatrixTable
    :param list bins: List of bins to filter by
    :return: Dict containing bin and remaning TP/FP count
    '''
    results = {}
    for bin in bins:
        print("Analysing bin " + str(bin))
        results[bin] = {}
        #filter by bin
        mt_tmp = mt.filter_rows(mt.info.bin <= bin)
        #genotype hard filters
        filter_condition = (
            (mt_tmp.GT.is_het() & (mt_tmp.HetAB < 0.3)) | 
            (mt_tmp.DP < 10) |
            (mt_tmp.GQ < 20)
        )
        mt_tmp = mt_tmp.annotate_entries(
            hard_filters = hl.if_else(filter_condition, 'Fail', 'Pass')
        )
        mt_tmp = mt_tmp.filter_entries(mt_tmp.hard_filters == 'Pass')
        #remove unused rows
        mt_tmp = hl.variant_qc(mt_tmp)
        mt_tmp = mt_tmp.filter_rows(mt_tmp.variant_qc.n_non_ref == 0, keep = False)
        #get counts
        counts = count_tp_fp(mt_tmp)
        results[bin]['TP'] = counts[0]
        results[bin]['FP'] = counts[1]

    return results


def filter_and_count(mt: hl.MatrixTable):
    '''
    Filter MT by various bins followed by genotype GQ and cauclate % of FP and TP remaining for each bin
    :param hl.MatrixTable mt: Input MatrixTable annotated with RF score, bin and TP/FP
    '''
    #split into SNPs and indels
    snp_mt = mt.filter_rows(hl.snp(mt.alleles[0], mt.alleles[1]))
    indel_mt = mt.filter_rows(hl.is_indel(mt.alleles[0], mt.alleles[1]))

    snp_total_tps, snp_total_fps = count_tp_fp(snp_mt)
    indel_total_tps, indel_total_fps = count_tp_fp(indel_mt)

    snp_bins = get_threshold_range(40)
    indel_bins = get_threshold_range(68)

    snp_results = filter_mt_count_tp_fp(snp_mt, snp_bins)
    indel_results = filter_mt_count_tp_fp(indel_mt, indel_bins)

    #calculate percentages and print results
    print("SNPs:")
    for bin in snp_results.keys():
        tp_pc = (snp_results[bin]['TP']/snp_total_tps)*100
        fp_pc = (snp_results[bin]['FP']/snp_total_fps)*100
        print("bin " + str(bin) + ": " + '{0:.3f}'.format(tp_pc) + " percent TPs remain, "
         + '{0:.3f}'.format(fp_pc) + " percent FPs remain" )

    print("Indels:")
    for bin in indel_results.keys():
        tp_pc = (indel_results[bin]['TP']/indel_total_tps)*100
        fp_pc = (indel_results[bin]['FP']/indel_total_fps)*100
        print("bin " + str(bin) + ": " + '{0:.3f}'.format(tp_pc) + " percent TPs remain, "
         + '{0:.3f}'.format(fp_pc) + " percent FPs remain" )



def main():
    #set up
    inputs = parse_config()
    mtdir = inputs['matrixtables_lustre_dir']
    rf_dir = inputs['var_qc_rf_dir']

    # initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    rf_htfile = rf_dir + "6617f838" + "/_gnomad_score_binning_tmp.ht"
    mtfile = mtdir + "mt_varqc_splitmulti.mt"

    mt = hl.read_matrix_table(mtfile)
    mt_annot = annotate_with_rf(mt, rf_htfile)
    filter_and_count(mt_annot)

if __name__ == '__main__':
    main()