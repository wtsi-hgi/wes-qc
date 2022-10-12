#compare different combinations of hard filters
import hail as hl
import pyspark
from wes_qc.utils.utils import parse_config


def annotate_with_rf(mt: hl.MatrixTable, rf_htfile: str) -> hl.MatrixTable:
    '''
    Annotate MatrixTable with TP, FP, rf_bin and rf_score
    :param hl.MatrixTable mt: Input MatrixTable
    :param str rf_htfile: Random forest ht file
    :return: hl.MatrixTable
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

    return mt


def annotate_cq(mt: hl.MatrixTable, cqfile: str) -> hl.MatrixTable:
    '''
    Annotate MatrixTable with consequence
    :param hl.MatrixTable mt: Input MatrixTable
    :param str cqfile: Most severe consequence annotation from VEP
    :return: hl.MatrixTable
    '''
    ht=hl.import_table(cqfile,types={'f0':'str','f1':'int32', 'f2':'str','f3':'str','f4':'str', 'f5':'str'}, no_header=True)
    ht=ht.annotate(chr=ht.f0)
    ht=ht.annotate(pos=ht.f1)
    ht=ht.annotate(rs=ht.f2)
    ht=ht.annotate(ref=ht.f3)
    ht=ht.annotate(alt=ht.f4)
    ht=ht.annotate(consequence=ht.f5)
    ht = ht.key_by(
    locus=hl.locus(ht.chr, ht.pos), alleles=[ht.ref,ht.alt])
    ht=ht.drop(ht.f0,ht.f1,ht.f2,ht.f3,ht.f4,ht.chr,ht.pos,ht.ref,ht.alt)
    ht = ht.key_by(ht.locus, ht.alleles)

    mt=mt.annotate_rows(consequence=ht[mt.row_key].consequence)

    return mt


def filter_and_count(mt: hl.MatrixTable, plot_dir: str, pedfile: str) -> dict:
    '''
    Filter MT by various bins followed by genotype GQ and cauclate % of FP and TP remaining for each bifn
    :param hl.MatrixTable mt: Input MatrixTable annotated with RF score, bin and TP/FP
    :param str plot_dir: directory for output files 
    :return: dict
    '''
    results = {'snv':{}, 'indel':{}}
    #split into SNPs and indels
    snp_mt = mt.filter_rows(hl.is_snp(mt.alleles[0], mt.alleles[1]))
    indel_mt = mt.filter_rows(hl.is_indel(mt.alleles[0], mt.alleles[1]))

    snp_total_tps, snp_total_fps = count_tp_fp(snp_mt)
    indel_total_tps, indel_total_fps = count_tp_fp(indel_mt)

    snp_bins = list(range(35,46))
    indel_bins = list(range(57,68))
    gq_vals = [20, 30]
    dp_vals = [10, 15, 20]
    ab_vals = [0.2, 0.3]

    for i_dp in dp_vals:
        dp = dp_vals[i_dp]
        dp_str = 'DP_' + str(dp)
        results['snv'][dp_str] = {}
        results['indel'][dp_str] = {}
        for i_gq in gq_vals:
            gq = gq_vals[i_gq]
            gq_str = 'GQ_' + str(gq)
            results['snv'][dp_str][gq_str] = {}
            results['indel'][dp_str][gq_str] = {}
            for i_ab in ab_vals:
                ab = ab_vals[i_ab]
                ab_str = 'AB_' + str(ab)
                snp_counts_per_bin = filter_mt_count_tp_fp_t_u(snp_mt, snp_bins, pedfile, dp, gq, ab, 'snv')
                indel_counts_per_bin = filter_mt_count_tp_fp_t_u(indel_mt, indel_bins, pedfile, dp, gq, ab, 'indel')
                results['snv'][dp_str][gq_str][ab_str] = snp_counts_per_bin
                results['indel'][dp_str][gq_str][ab_str] = indel_counts_per_bin
                print(results)
                exit(0)


def filter_mt_count_tp_fp_t_u(mt: hl.MatrixTable, bins: list, pedfile: str, dp: int, gq: int, ab: float, var_type: str):
    '''
    Filter mt by each rf bin in a list, then genotype hard filters and count remaining TP and P variants
    :param hl.MatrixTable mt: Input MatrixTable
    :param list bins: List of bins to filter by
    :param str pedfile: pedfile path
    :param int dp: DP threshold
    ;param int gq; GQ threshold
    :param float ab: allele balance threshold
    :param  str var_type: variant type (snv/indel)
    :return: Dict containing bin and remaning TP/FP count
    '''
    results = {}
    pedigree = hl.Pedigree.read(pedfile)
    #list of samples in trios
    trio_sample_ht = hl.import_fam(pedfile)
    sample_list = trio_sample_ht.id.collect() + trio_sample_ht.pat_id.collect() + trio_sample_ht.mat_id.collect()
    
    for bin in bins:
        print("Analysing bin " + str(bin))
        results[bin] = {}
        #filter by bin
        mt_tmp = mt.filter_rows(mt.info.rf_bin <= bin)
        #genotype hard filters
        filter_condition = (
            (mt_tmp.GT.is_het() & (mt_tmp.HetAB < ab)) | 
            (mt_tmp.DP < dp) |
            (mt_tmp.GQ < gq)
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
        #get transmitted/untransmitted
        if var_type == 'snv':
            ratio = get_trans_untrans(mt_tmp, pedigree, sample_list)
            results[bin]['t_u_ratio'] = ratio
        else:
            results[bin]['t_u_ratio'] = 0

    return results


def get_trans_untrans(mt: hl.MatrixTable, pedigree: hl.Pedfile, sample_list: list) -> float:
    '''
    get transmitted/untransmitted ratio
    :param hl.MatrixTable mt: matrixtable
    :param hl.Pedigree pedigree: Hail Pedigree
    :param list sample_list: List of samples in trios
    :return float:
    '''
    #filter to synonymous
    mt_syn = mt.filter_rows(mt.consequence == 'synonymous_variant')
        #restrict to samples in trios, annotate with AC and filter to trio AC == 1 or 2
    mt2 = mt_syn.filter_cols(hl.set(sample_list).contains(mt_syn.s))
    mt2 = hl.variant_qc(mt2, name='varqc_trios')
        #split to potentially transitted/untransmitted
    untrans_mt = mt2.filter_rows(mt2.varqc_trios.AC[1] == 1)
    trans_mt = mt2.filter_rows(mt2.varqc_trios.AC[1] == 2)
        #run tdt function for potential trans and untrans
    tdt_ht_trans = hl.transmission_disequilibrium_test(trans_mt, pedigree)
    tdt_ht_untrans = hl.transmission_disequilibrium_test(untrans_mt, pedigree)
    trans_sing = tdt_ht_trans.filter((tdt_ht_trans.t == 1) & (tdt_ht_trans.u == 0))
    trans = trans_sing.count()
    untrans_sing = tdt_ht_untrans.filter((tdt_ht_untrans.t == 0) & (tdt_ht_untrans.u == 1))
    untrans = untrans_sing.count()
    if untrans > 0:
        ratio = trans/untrans
    else:
        ratio = 0

    return ratio


def count_tp_fp(mt: hl.MatrixTable) -> tuple:
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


def main():
    # set up
    inputs = parse_config()
    mtdir = inputs['matrixtables_lustre_dir']
    rf_dir = inputs['var_qc_rf_dir']
    resourcedir = inputs['resource_dir']
    plot_dir = inputs['plots_dir_local']

    # initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    rf_htfile = rf_dir + "6617f838" + "/_gnomad_score_binning_tmp.ht"
    mtfile = mtdir + "mt_varqc_splitmulti.mt"
    cqfile = resourcedir + "all_consequences.txt"
    pedfile = resourcedir + "trios.ped"

    mt = hl.read_matrix_table(mtfile)
    mt_annot = annotate_with_rf(mt, rf_htfile)
    mt_annot = annotate_cq(mt_annot, cqfile)
    results = filter_and_count(mt_annot, plot_dir, pedfile)


if __name__ == '__main__':
    main()