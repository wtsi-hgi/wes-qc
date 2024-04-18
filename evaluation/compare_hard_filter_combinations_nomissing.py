#compare different combinations of hard filters
import hail as hl
import pyspark
import datetime
import argparse
from wes_qc.utils.utils import parse_config


def get_options():
    '''
    Get options from the command line
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("--runhash", help="RF run hash")
    parser.add_argument("--generate-giab", help="Re-generate annotated GIAB table", action="store_true")
    args = parser.parse_args()
    if not args.runhash:
        print("--runhash must be specified")
        exit(1)

    return args


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


def prepare_giab_ht(giab_vcf: str, giab_cqfile: str, mtdir: str) -> hl.Table:
    '''
    Get GIAB ht from vcf file
    :param str giab_vcf: path of input VCF file
    :param str giab_cqfile: path of vep annotation
    :param str mtdir: MatrixTable directory
    :return: hl.Table
    '''
    mt = hl.import_vcf(giab_vcf, force_bgz = True, reference_genome='GRCh38')
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(mt.variant_qc.n_non_ref > 0)

    ht=hl.import_table(giab_cqfile,types={'f0':'str','f1':'int32', 'f2':'str','f3':'str','f4':'str', 'f5':'str', 'f6':'str', 'f7': 'str'}, no_header=True)
    ht=ht.annotate(chr=ht.f0)
    ht=ht.annotate(pos=ht.f1)
    ht=ht.annotate(rs=ht.f2)
    ht=ht.annotate(ref=ht.f3)
    ht=ht.annotate(alt=ht.f4)
    ht=ht.annotate(consequence=ht.f5)
    ht=ht.annotate(impacte=ht.f6)
    ht=ht.annotate(misc=ht.f7)
    ht = ht.key_by(
    locus=hl.locus(ht.chr, ht.pos), alleles=[ht.ref,ht.alt])
    ht=ht.drop(ht.f0,ht.f1,ht.f2,ht.f3,ht.f4,ht.chr,ht.pos,ht.ref,ht.alt)
    ht = ht.key_by(ht.locus, ht.alleles)

    mt=mt.annotate_rows(consequence=ht[mt.row_key].consequence)
    giab_vars = mt.rows()
    tmphtg = mtdir + "tmphtgx.ht"
    giab_vars = giab_vars.checkpoint(tmphtg, overwrite = True)

    return giab_vars


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


def filter_and_count(mt_tp: hl.MatrixTable, mt_fp: hl.MatrixTable, mt_syn: hl.MatrixTable, mt_prec_recall: hl.Table, ht_giab: hl.Table, plot_dir: str, pedfile: str, mtdir: str) -> dict:
    '''
    Filter MT by various bins followed by genotype GQ and cauclate % of FP and TP remaining for each bifn
    :param hl.MatrixTable mt_tp: Input TP MatrixTable
    :param hl.MatrixTable mt_fp: Input FP MatrixTable
    :param hl.MatrixTable mt_syn: Input synonymous MatrixTable
    :param hl.MatrixTable mt_prec_recall: mt from GIAB sample for precision/recall analysys
    :param hl.Table ht_giab: GIAB variants
    :param str plot_dir: directory for output files 
    :param str mtdir: matrixtable directory
    :return: dict
    '''
    results = {'snv':{}, 'indel':{}}
    #split mts into SNPs and indels
    snp_mt_tp = mt_tp.filter_rows(hl.is_snp(mt_tp.alleles[0], mt_tp.alleles[1]))
    indel_mt_tp = mt_tp.filter_rows(hl.is_indel(mt_tp.alleles[0], mt_tp.alleles[1]))
    snp_mt_fp = mt_fp.filter_rows(hl.is_snp(mt_fp.alleles[0], mt_fp.alleles[1]))
    indel_mt_fp = mt_fp.filter_rows(hl.is_indel(mt_fp.alleles[0], mt_fp.alleles[1]))

    snp_total_tps, snp_total_fps = count_tp_fp(snp_mt_tp, snp_mt_fp)
    indel_total_tps, indel_total_fps = count_tp_fp(indel_mt_tp, indel_mt_fp)

    results['snv_total_tp'] = snp_total_tps
    results['snv_total_fp'] = snp_total_fps
    results['indel_total_tp'] = indel_total_tps
    results['indel_total_fp'] = indel_total_fps

    snp_bins = list(range(35,46))
    indel_bins = list(range(58,69))
    gq_vals = [10, 15, 20]
    # dp_vals = [4, 5, 6, 10]
    # ab_vals = [0.2, 0.25, 0.3]
    dp_vals = [5, 10]
    ab_vals = [0.2, 0.3]

    for bin in snp_bins:
        bin_str = "bin_" + str(bin)
        print(f"= Processing SNP bin {bin}")
        mt_tp_tmp = snp_mt_tp.filter_rows(snp_mt_tp.info.rf_bin <= bin)
        tmpmtb1 = mtdir + "tmp1bx.mt"
        mt_tp_tmp = mt_tp_tmp.checkpoint(tmpmtb1, overwrite = True)
        mt_fp_tmp = snp_mt_fp.filter_rows(snp_mt_fp.info.rf_bin <= bin)
        tmpmtb2 = mtdir + "tmp2bx.mt"
        mt_fp_tmp = mt_fp_tmp.checkpoint(tmpmtb2, overwrite = True)
        mt_syn_tmp = mt_syn.filter_rows(mt_syn.info.rf_bin <= bin)
        tmpmtb3 = mtdir + "tmp3bx.mt"
        mt_syn_tmp = mt_syn_tmp.checkpoint(tmpmtb3, overwrite = True)
        mt_prec_recall_tmp = mt_prec_recall.filter_rows(mt_prec_recall.info.rf_bin <= bin)
        tmphtb4 = mtdir + "tmp4bx.ht"
        mt_prec_recall_tmp = mt_prec_recall_tmp.checkpoint(tmphtb4, overwrite = True)

        for dp in dp_vals:
            dp_str = 'DP_' + str(dp)
            for gq in gq_vals:
                gq_str = 'GQ_' + str(gq)
                for ab in ab_vals:
                    ab_str = 'AB_' + str(ab)
                    print("-- Filter combination: " + dp_str + " " + gq_str + " " + ab_str)
                    filter_name = ("_").join([bin_str, dp_str, gq_str, ab_str])
                    snp_counts = filter_mt_count_tp_fp_t_u(mt_tp_tmp, mt_fp_tmp, mt_syn_tmp, mt_prec_recall_tmp, ht_giab, pedfile, dp, gq, ab, 'snv', mtdir)
                    results['snv'][filter_name] = snp_counts

    for bin in indel_bins:
        bin_str = "bin_" + str(bin)
        print(f"= Processing InDel bin :{bin}")
        mt_tp_tmp = indel_mt_tp.filter_rows(indel_mt_tp.info.rf_bin <= bin)
        tmpmtb1 = mtdir + "tmp1bx.mt"
        mt_tp_tmp = mt_tp_tmp.checkpoint(tmpmtb1, overwrite = True)
        mt_fp_tmp = indel_mt_fp.filter_rows(indel_mt_fp.info.rf_bin <= bin)
        tmpmtb2 = mtdir + "tmp2bx.mt"
        mt_fp_tmp = mt_fp_tmp.checkpoint(tmpmtb2, overwrite = True)
        mt_prec_recall_tmp = mt_prec_recall.filter_rows(mt_prec_recall.info.rf_bin <= bin)
        tmphtb4 = mtdir + "tmp4bx.ht"
        mt_prec_recall_tmp = mt_prec_recall_tmp.checkpoint(tmphtb4, overwrite = True)

        for dp in dp_vals:
            dp_str = 'DP_' + str(dp)
            for gq in gq_vals:
                gq_str = 'GQ_' + str(gq)
                for ab in ab_vals:
                    ab_str = 'AB_' + str(ab)
                    print("-- Filter combination: " + dp_str + " " + gq_str + " " + ab_str)
                    filter_name = ("_").join([bin_str, dp_str, gq_str, ab_str])
                    indel_counts = filter_mt_count_tp_fp_t_u(mt_tp_tmp, mt_fp_tmp, mt_syn, mt_prec_recall_tmp, ht_giab, pedfile, dp, gq, ab, 'indel', mtdir)
                    results['indel'][filter_name] = indel_counts



    return results


def filter_mt_count_tp_fp_t_u(mt_tp: hl.MatrixTable, mt_fp: hl.MatrixTable, mt_syn: hl.MatrixTable, mt_prec_recall: hl.Table, ht_giab: hl.Table, pedfile: str, dp: int, gq: int, ab: float, var_type: str, mtdir: str):
    '''
    Filter mt by each rf bin in a list, then genotype hard filters and count remaining TP and P variants
    :param hl.MatrixTable mt_tp: Input TP MatrixTable
    :param hl.MatrixTable mt_fp: Input FP MatrixTable
    :param hl.MatrixTable mt_syn: Input synonymous MatrixTable
    :param hl.MatrixTable mt_prec_recall: mt from GIAB sample for precision/recall
    :param hl.Table ht_giab: GIAB variants
    :param str pedfile: pedfile path
    :param int dp: DP threshold
    ;param int gq; GQ threshold
    :param float ab: allele balance threshold
    :param  str var_type: variant type (snv/indel)
    :param str mtdir: matrixtable directory
    :return: Dict containing bin and remaning TP/FP count
    '''
    results = {}
    pedigree = hl.Pedigree.read(pedfile)
    #list of samples in trios
    trio_sample_ht = hl.import_fam(pedfile)
    sample_list = trio_sample_ht.id.collect() + trio_sample_ht.pat_id.collect() + trio_sample_ht.mat_id.collect()

    #genotype hard filters - should put this in a different method
    now = datetime.datetime.now()
    print(now.time())
    filter_condition = (
        (mt_tp.GT.is_het() & (mt_tp.HetAB < ab)) | 
        (mt_tp.DP < dp) |
        (mt_tp.GQ < gq)
    )
    mt_tp_tmp = mt_tp.annotate_entries(
        hard_filters = hl.if_else(filter_condition, 'Fail', 'Pass')
    )
    mt_tp_tmp = mt_tp_tmp.filter_entries(mt_tp_tmp.hard_filters == 'Pass')
        #remove unused rows
    mt_tp_tmp = hl.variant_qc(mt_tp_tmp)
    mt_tp_tmp = mt_tp_tmp.filter_rows(mt_tp_tmp.variant_qc.n_non_ref == 0, keep = False)

    filter_condition = (
        (mt_fp.GT.is_het() & (mt_fp.HetAB < ab)) | 
        (mt_fp.DP < dp) |
        (mt_fp.GQ < gq)
    )
    mt_fp_tmp = mt_fp.annotate_entries(
        hard_filters = hl.if_else(filter_condition, 'Fail', 'Pass')
    )
    mt_fp_tmp = mt_fp_tmp.filter_entries(mt_fp_tmp.hard_filters == 'Pass')
        #remove unused rows
    mt_fp_tmp = hl.variant_qc(mt_fp_tmp)
    mt_fp_tmp = mt_fp_tmp.filter_rows(mt_fp_tmp.variant_qc.n_non_ref == 0, keep = False)

    filter_condition = (
        (mt_syn.GT.is_het() & (mt_syn.HetAB < ab)) | 
        (mt_syn.DP < dp) |
        (mt_syn.GQ < gq)
    )
    mt_syn_tmp = mt_syn.annotate_entries(
        hard_filters = hl.if_else(filter_condition, 'Fail', 'Pass')
    )
    mt_syn_tmp = mt_syn_tmp.filter_entries(mt_syn_tmp.hard_filters == 'Pass')
        #remove unused rows
    mt_syn_tmp = hl.variant_qc(mt_syn_tmp)
    mt_syn_tmp = mt_syn_tmp.filter_rows(mt_syn_tmp.variant_qc.n_non_ref == 0, keep = False)

    filter_condition = (
        (mt_prec_recall.GT.is_het() & (mt_prec_recall.HetAB < ab)) | 
        (mt_prec_recall.DP < dp) |
        (mt_prec_recall.GQ < gq)
    )
    mt_prec_recall_tmp = mt_prec_recall.annotate_entries(
        hard_filters = hl.if_else(filter_condition, 'Fail', 'Pass')
    )
    mt_prec_recall_tmp = mt_prec_recall_tmp.filter_entries(mt_prec_recall_tmp.hard_filters == 'Pass')
        #remove unused rows
    mt_prec_recall_tmp = hl.variant_qc(mt_prec_recall_tmp)
    mt_prec_recall_tmp = mt_prec_recall_tmp.filter_rows(mt_prec_recall_tmp.variant_qc.n_non_ref == 0, keep = False)
    ht_prec_recall_tmp = mt_prec_recall_tmp.rows()

    # tmpmt2 = mtdir + "tmp2.mt"
    # mt_tmp = mt_tmp.checkpoint(tmpmt2, overwrite = True)
    counts = count_tp_fp(mt_tp_tmp, mt_fp_tmp)
    results['TP'] = counts[0]
    results['FP'] = counts[1]
    if var_type == 'snv':
        ratio = get_trans_untrans(mt_syn_tmp, pedigree, sample_list, mtdir)
        prec, recall = get_prec_recall(ht_prec_recall_tmp, ht_giab, 'snv', mtdir)
        results['t_u_ratio'] = ratio
        results['prec'] = prec
        results['recall'] = recall
    else:
        prec, recall, prec_frameshift, recall_frameshift, prec_inframe, recall_inframe = get_prec_recall(ht_prec_recall_tmp, ht_giab, 'indel', mtdir)
        results['prec'] = prec
        results['recall'] = recall
        results['prec_inframe'] = prec_inframe
        results['recall_inframe'] = recall_inframe
        results['prec_frameshift'] = prec_frameshift
        results['recall_frameshift'] = recall_frameshift

    return results


def get_prec_recall(ht_prec_recall: hl.Table, ht_giab: hl.Table, var_type: str, mtdir: str) -> tuple:
    '''
    Get precison/recall vs GIAB
    :param hl.Table ht_prec_recall: hail Table of variants in ALSPAC GIAB sample
    :param hl.Table ht_giab: hail Table GIAB variants
    :param str var_type: Variant type
    :param str mtdir: matrixtable directory
    :return: tuple
    '''

    if var_type == 'snv':
        giab_snvs = ht_giab.filter(hl.is_snp(ht_giab.alleles[0], ht_giab.alleles[1]))
        alspac_snvs = ht_prec_recall.filter(hl.is_snp(ht_prec_recall.alleles[0], ht_prec_recall.alleles[1]))
        tmpht1 = mtdir + "tmppr1x.ht"
        giab_snvs = giab_snvs.checkpoint(tmpht1, overwrite = True)
        tmpht2 = mtdir + "tmppr2x.ht"
        alspac_snvs = alspac_snvs.checkpoint(tmpht2, overwrite = True)
        p, r = calculate_precision_recall(giab_snvs, alspac_snvs)
        return p, r

    elif var_type == 'indel':
        giab_indels = ht_giab.filter(hl.is_indel(ht_giab.alleles[0], ht_giab.alleles[1]))
        alspac_indels = ht_prec_recall.filter(hl.is_indel(ht_prec_recall.alleles[0], ht_prec_recall.alleles[1]))
        tmpht1 = mtdir + "tmppr1x.ht"
        giab_indels = giab_indels.checkpoint(tmpht1, overwrite = True)
        tmpht2 = mtdir + "tmppr2x.ht"
        alspac_indels = alspac_indels.checkpoint(tmpht2, overwrite = True)
        p, r = calculate_precision_recall(giab_indels, alspac_indels)

        giab_frameshift = giab_indels.filter(giab_indels.consequence == 'frameshift_variant')
        alspac_frameshift = alspac_indels.filter(alspac_indels.consequence == 'frameshift_variant')
        tmpht3 = mtdir + "tmppr3x.ht"
        giab_frameshift = giab_frameshift.checkpoint(tmpht3, overwrite = True)
        tmpht4 = mtdir + "tmppr4x.ht"
        alspac_frameshift = alspac_frameshift.checkpoint(tmpht4, overwrite = True)
        p_f, r_f = calculate_precision_recall(giab_frameshift, alspac_frameshift)

        inframe_cqs = ['inframe_deletion', 'inframe_insertion']
        giab_in_frame = giab_indels.filter(hl.literal(inframe_cqs).contains(giab_indels.consequence))
        alspac_in_frame = alspac_indels.filter(hl.literal(inframe_cqs).contains(alspac_indels.consequence))
        tmpht5 = mtdir + "tmppr5x.ht"
        giab_in_frame = giab_in_frame.checkpoint(tmpht5, overwrite = True)
        tmpht6 = mtdir + "tmppr6x.ht"
        alspac_in_frame = alspac_in_frame.checkpoint(tmpht6, overwrite = True)
        p_if, r_if = calculate_precision_recall(giab_in_frame, alspac_in_frame)

        return p, r, p_f, r_f, p_if, r_if


def calculate_precision_recall(ht_control: hl.Table, ht_test: hl.Table) -> tuple:
    '''
    Calculate orecision recall
    :param hl.Table ht_control: Control set
    :paran hl.Table ht_test: Test set
    :return tuple:
    '''
    print("get intersects")
    vars_in_both = ht_control.semi_join(ht_test)
    control_only = ht_control.anti_join(ht_test)
    test_only = ht_test.anti_join(ht_control)
    print("count_vars")
    tp = vars_in_both.count()
    fn = control_only.count()
    fp = test_only.count()

    precision = tp / (tp + fp)
    recall = tp / (tp + fn)

    return precision, recall


def get_trans_untrans(mt: hl.MatrixTable, pedigree: hl.Pedigree, sample_list: list, mtdir: str) -> float:
    '''
    get transmitted/untransmitted ratio
    :param hl.MatrixTable mt: matrixtable
    :param hl.Pedigree pedigree: Hail Pedigree
    :param list sample_list: List of samples in trios
    :param str mtdir: matrixtable directory
    :return float:
    '''
    #filter to synonymous
    mt_syn = mt.filter_rows(mt.consequence == 'synonymous_variant')
        #restrict to samples in trios, annotate with AC and filter to trio AC == 1 or 2
    mt2 = mt_syn.filter_cols(hl.set(sample_list).contains(mt_syn.s))
    mt2 = hl.variant_qc(mt2, name='varqc_trios')
    tmpmt3 = mtdir + "tmp3x.mt"
    mt2 = mt2.checkpoint(tmpmt3, overwrite = True)
    #split to potentially transitted/untransmitted
    untrans_mt = mt2.filter_rows(mt2.varqc_trios.AC[1] == 1)
    tmpmt4 = mtdir + "tmp4x.mt"
    untrans_mt = untrans_mt.checkpoint(tmpmt4, overwrite = True)
    trans_mt = mt2.filter_rows(mt2.varqc_trios.AC[1] == 2)
    tmpmt5 = mtdir + "tmp5x.mt"
    trans_mt = trans_mt.checkpoint(tmpmt5, overwrite = True)
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


def count_tp_fp(mt_tp: hl.MatrixTable, mt_fp: hl.MatrixTable) -> tuple:
    '''
    Count total TPs and FPs in a pair of MatrixTables
    :param hl.MatrixTable mt_tp: Input TP mt
    :param hl.MatrixTable mt_fp: Input FP mt
    :return: tuple of TP and FP counts
    '''

    tp_count = mt_tp.count_rows()
    fp_count = mt_fp.count_rows()

    return tp_count, fp_count


def filter_mts(mt: hl.MatrixTable, mtdir: str) -> tuple:
    '''
    Split matrixtable and return tables with just TP, just FP, just synonymous
    :param hl.MatrixTable mt: Input mtfile
    :param st mtdir: matrixtable directory
    :return: tuple of 4 hl.MatrixTable objects
    '''
    mt_true = mt.filter_rows(mt.TP == True)#TP variants
    mt_false = mt.filter_rows(mt.FP == True)#FP variants
    mt_syn = mt.filter_rows(mt.consequence == 'synonymous_variant')#synonymous for transmitted/unstransmitted
    sample = 'EGAN00003332049'#GIAB12878/HG001
    mt_prec_recall = mt.filter_cols(mt.s == sample)#GIAB sample for precision/recall
    mt_prec_recall = mt_prec_recall.filter_rows(mt_prec_recall.locus.in_autosome())
    mt_prec_recall = hl.variant_qc(mt_prec_recall)
    mt_prec_recall = mt_prec_recall.filter_rows(mt_prec_recall.variant_qc.n_non_ref > 0)
    tmpmtt = mtdir + "tmptx.mt"
    tmpmtf = mtdir + "tmpfx.mt"
    tmpmts = mtdir + "tmpsx.mt"
    tmpmtpr = mtdir + "tmpprx.mt"
    mt_true = mt_true.checkpoint(tmpmtt, overwrite = True)
    mt_false = mt_false.checkpoint(tmpmtf, overwrite = True)
    mt_syn = mt_syn.checkpoint(tmpmts, overwrite = True)
    mt_prec_recall = mt_prec_recall.checkpoint(tmpmtpr, overwrite = True)

    return mt_true, mt_false, mt_syn, mt_prec_recall


def print_results(results: dict, outfile: str, vartype: str):
    '''
    Print results dict to a file
    :param dict results: results dict
    :param str outfile: output file path
    :param str vartype: variant type (snv or indel)
    '''
    header = ['filter', 'TP', 'FP']
    if vartype == 'snv':
        header = header + (['t_u_ratio', 'precision', 'recall'])
    elif vartype == 'indel':
        header = header + (['precision', 'recall', 'precision_frameshift', 'recall_frameshift', 'precision_inframe', 'recall_inframe'])
    
    with open(outfile, 'w') as o:
        o.write(("\t").join(header))
        o.write("\n")

        for var_f in results[vartype].keys():
            if vartype == 'snv':
                tp = str((results[vartype][var_f]['TP']/results['snv_total_tp']) * 100)
                fp = str((results[vartype][var_f]['FP']/results['snv_total_fp']) * 100)
            elif vartype == 'indel':
                tp = str((results[vartype][var_f]['TP']/results['indel_total_tp']) * 100)
                fp = str((results[vartype][var_f]['FP']/results['indel_total_fp']) * 100)
            outline = [var_f, tp, fp]
            if vartype == 'snv':
                tu = str(results[vartype][var_f]['t_u_ratio'])
                p = str(results[vartype][var_f]['prec'])
                r = str(results[vartype][var_f]['recall'])
                outline = outline + [tu, p, r]
            elif vartype == 'indel':
                p = str(results[vartype][var_f]['prec'])
                r = str(results[vartype][var_f]['recall']) 
                p_f = str(results[vartype][var_f]['prec_frameshift'])
                r_f = str(results[vartype][var_f]['recall_frameshift']) 
                p_if = str(results[vartype][var_f]['prec_inframe'])
                r_if = str(results[vartype][var_f]['recall_inframe']) 
                outline = outline + [p, r, p_f, r_f, p_if, r_if]

            o.write(("\t").join(outline))
            o.write("\n")


def main():
    # set up
    args = get_options()
    inputs = parse_config()
    mtdir = inputs['matrixtables_lustre_dir']
    rf_dir = inputs['var_qc_rf_dir']
    resourcedir = inputs['resource_dir']
    plot_dir = inputs['plots_dir_local']

    # initialise hail
    tmp_dir = inputs["tmp_dir"]
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    giab_ht_path = mtdir + "giab_annotated.ht"
    if args.generate_giab:
        print("Generating GIAB matrix")
        giab_vcf = resourcedir + "HG001_GRCh38_benchmark.interval.illumina.vcf.gz"
        giab_cqfile = resourcedir + "all.interval.illumina.vep.info.txt"
        giab_ht = prepare_giab_ht(giab_vcf, giab_cqfile, mtdir)
        giab_ht.write(giab_ht_path)
    else:
        giab_ht = hl.read_table(giab_ht_path)

    rf_htfile = rf_dir + args.runhash + "/_gnomad_score_binning_tmp.ht"
    mtfile = mtdir + "mt_varqc_splitmulti.mt"
    cqfile = resourcedir + "all_consequences.txt"
    pedfile = resourcedir + "trios.ped"


    mt = hl.read_matrix_table(mtfile)
    print(f"Loaded samples matrix {mt.count_cols()} samples, {mt.count_rows()} variants")

    print("Annotating matrixtable")
    mt_annot = annotate_with_rf(mt, rf_htfile)
    mt_annot = annotate_cq(mt_annot, cqfile)
    mt_tp, mt_fp, mt_syn, mt_prec_recall = filter_mts(mt_annot, mtdir)

    print("=== Starting to search for hard filter combinations ===")
    results = filter_and_count(mt_tp, mt_fp, mt_syn, mt_prec_recall, giab_ht, plot_dir, pedfile, mtdir)

    outfile_snv = plot_dir + "/" + args.runhash + "_genotype_hard_filter_comparison_snv_fewer_combs_v2.txt"
    outfile_indel = plot_dir + "/" + args.runhash + "_genotype_hard_filter_comparison_indel_fewer_combs_v2.txt"
    print_results(results, outfile_snv, 'snv')
    print_results(results, outfile_indel, 'indel')


if __name__ == '__main__':
    main()