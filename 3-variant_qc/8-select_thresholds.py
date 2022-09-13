# select thresholds for random forest - for a range of thresholds for SNPs and indels 
# calculate how many TPs are retained and how many FPs are included
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
    parser.add_argument("--snv", help = "Proposed SNV threshold")
    parser.add_argument("--indel", help = "Proposed Indel threshold")
    args = parser.parse_args()
    if not args.runhash and args.snv and args.indel:
        print("--runhash, --snv and --indel must be specified")
        exit(1)
    args.snv = int(args.snv)
    args.indel = int(args.indel)

    return args


def analyse_thresholds(htfile: str, snv_threshold: int, indel_threshold: int):
    '''
    For suggested thresholds calculate what proportion of TP variants are included and what propeortion of FPs
    :param str htfile: ht file annotated with rank, bin, tp and fp
    :param int snv_threshold: proposed threshold for SNVs
    :param int indel_threshold: proposed threshold for indels
    '''
    #for each threshold also look at 3 bins on either side (unless <1 or >101)
    snv_t = get_threshold_range(snv_threshold)
    indel_t = get_threshold_range(indel_threshold)
    
    ht = hl.read_table(htfile)
    #keep only those with rank_id = rank
    ht = ht.filter(ht.rank_id == 'rank')
    #split by snv and indel
    snv_ht = ht.filter(ht.allele_type == 'snv')
    indel_ht = ht.filter(ht.allele_type == 'snv', keep = False)

    #calculate % TP kept and % FP kept at each threshold
    snv_result = get_vars_kept(snv_ht, snv_t)
    indel_result = get_vars_kept(indel_ht, indel_t)
    #print results
    print("SNV percentage of variants kept at each threshold:")
    print_results(snv_result)
    print("Indel percentage of variants kept at each threshold:")
    print_results(indel_result)


def print_results(result: dict):
    '''
    Print percentages of variants kept at each threshold
    :param dict result: Results to print
    '''
    for t in result.keys():
        print("bin " + str(t) + ": " + '{0:.3f}'.format(result[t]['tp_pc']) + " percent TPs remain, "
         + '{0:.3f}'.format(result[t]['fp_pc']) + " percent FPs remain" )


def get_vars_kept(ht: hl.Table, t_list: list):
    '''
    find proportions fo TP and FPs kept in ht for each threshold in t_list
    :param hl.Table ht: input Hail Table
    :param list t_list: List of thresholds
    '''
    results = {}
    ht_TP = ht.filter(ht.tp == True)
    ht_FP = ht.filter(ht.fp == True)
    tp_total = ht_TP.count()
    fp_total = ht_FP.count()
    for t in t_list:
        results[t] = {}
        ht_tp_pass = ht_TP.filter(ht_TP.bin <= t)
        tp_pass = ht_tp_pass.count()
        results[t]['tp_pc'] = 100*(tp_pass/tp_total)
        ht_fp_pass = ht_TP.filter(ht_FP.bin <= t)
        fp_pass = ht_fp_pass.count()
        results[t]['fp_pc'] = 100*(fp_pass/fp_total)

    return results
    

def get_threshold_range(threshold: int):
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


def main():
    # set up
    args = get_options()
    inputs = parse_config()
    rf_dir = inputs['var_qc_rf_dir']

    # initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    htfile = rf_dir + args.runhash + "/_gnomad_score_binning_tmp.ht"
    analyse_thresholds(htfile, args.snv, args.indel)


if __name__ == '__main__':
    main()