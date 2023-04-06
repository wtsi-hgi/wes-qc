import json
import os.path

import hail as hl
import pyspark
import datetime
from utils.utils import parse_config
from typing import Tuple

from evaluation.roh_compare import count_hets_in_rohs, import_roh, write_out
from evaluation.compare_hard_filter_combinations import annotate_with_rf
from evaluation.compare_hard_filter_combinations import annotate_cq
from evaluation.compare_hard_filter_combinations import count_tp_fp
from evaluation.compare_hard_filter_combinations import get_trans_untrans

rf_hash = 'ce85819e'
snp_bins = [78, 81]
indel_bins = [49, 52, 55]
gq_vals = [10, 15, 20]
dp_vals = [5, 10]
ab_vals = [0.2, 0.3]
roh_path = 'file:///lustre/scratch123/projects/gnh_industry/Genes_and_Health_2023_02_44k/GH_44k_autosome_maf0.01_geno0.01_hwe1e-6_ROH_CALLING_OUT.hail-ready.hom'

snp_label = 'snp'
indel_label = 'indel'


def clean_mt(mt: hl.MatrixTable) -> hl.MatrixTable:
    mt = mt.select_entries(mt.GT, mt.HetAB, mt.DP, mt.GQ)
    mt = mt.drop(mt.assigned_pop, *mt.row_value)
    mt = mt.annotate_rows(info=hl.Struct())
    mt = mt.annotate_rows(
        type=hl.case()
               .when(hl.is_snp(mt.alleles[0], mt.alleles[1]), snp_label)
               .when(hl.is_indel(mt.alleles[0], mt.alleles[1]), indel_label)
               .default('other')
    )
    return mt


def filter_mts(mt: hl.MatrixTable, roh_path: str, mtdir: str) -> Tuple[str, str, str, str]:
    """
    Split matrixtable and return tables with just TP, just FP, just synonymous
    :param hl.MatrixTable mt: Input mtfile
    :param st mtdir: matrixtable directory
    :return: tuple of 4 hl.MatrixTable objects
    """
    mt_true = mt.filter_rows(mt.TP == True)  # TP variants
    mt_false = mt.filter_rows(mt.FP == True)  # FP variants
    mt_syn = mt.filter_rows(mt.consequence == 'synonymous_variant')  # synonymous for transmitted/unstransmitted

    tmpmtt = os.path.join(mtdir, "tp.mt")
    tmpmtf = os.path.join(mtdir, "fp.mt")
    tmpmts = os.path.join(mtdir, "syn.mt")
    tmpmtr = os.path.join(mtdir, "roh.mt")

    mt_true.write(tmpmtt, overwrite=True)
    mt_false.write(tmpmtf, overwrite=True)
    mt_syn.write(tmpmts, overwrite=True)

    roh = import_roh(roh_path)
    samples = roh.aggregate(hl.agg.collect_as_set(roh.ID))
    mt_roh = mt.filter_cols(hl.set(samples).contains(mt.s))
    mt_roh.write(tmpmtr, overwrite=True)

    return tmpmtt, tmpmtf, tmpmts, tmpmtr


def count_rows(mt: hl.MatrixTable, mutation_type: str) -> int:
    c = mt.filter_rows(mt.type == mutation_type).count_rows()
    return c


def apply_hard_filters(mt: hl.MatrixTable, dp: int, gq: int, ab: float) -> hl.MatrixTable:
    filter_condition = (
            (mt.GT.is_het() & (mt.HetAB < ab)) |
            (mt.DP < dp) |
            (mt.GQ < gq)
    )
    mt_tmp = mt.annotate_entries(
        hard_filters=hl.if_else(filter_condition, 'Fail', 'Pass')
    )
    mt_tmp = mt_tmp.filter_entries(mt_tmp.hard_filters == 'Pass')
    return mt_tmp


def filter_mt_count_tp_fp_t_u(mt_tp: hl.MatrixTable, mt_fp: hl.MatrixTable, mt_syn: hl.MatrixTable, pedfile: str, dp: int, gq: int, ab: float, var_type: str, mtdir: str):
    '''
    Filter mt by each rf bin in a list, then genotype hard filters and count remaining TP and P variants
    :param hl.MatrixTable mt_tp: Input TP MatrixTable
    :param hl.MatrixTable mt_fp: Input FP MatrixTable
    :param hl.MatrixTable mt_syn: Input synonymous MatrixTable
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
    mt_tp_tmp = apply_hard_filters(mt_tp, dp=dp, gq=gq, ab=ab)
        #remove unused rows
    mt_tp_tmp = hl.variant_qc(mt_tp_tmp)
    mt_tp_tmp = mt_tp_tmp.filter_rows(mt_tp_tmp.variant_qc.n_non_ref == 0, keep = False)

    mt_fp_tmp = apply_hard_filters(mt_fp, dp=dp, gq=gq, ab=ab)
        #remove unused rows
    mt_fp_tmp = hl.variant_qc(mt_fp_tmp)
    mt_fp_tmp = mt_fp_tmp.filter_rows(mt_fp_tmp.variant_qc.n_non_ref == 0, keep = False)

    mt_syn_tmp = apply_hard_filters(mt_syn, dp=dp, gq=gq, ab=ab)
        #remove unused rows
    mt_syn_tmp = hl.variant_qc(mt_syn_tmp)
    mt_syn_tmp = mt_syn_tmp.filter_rows(mt_syn_tmp.variant_qc.n_non_ref == 0, keep = False)

    # tmpmt2 = mtdir + "tmp2.mt"
    # mt_tmp = mt_tmp.checkpoint(tmpmt2, overwrite = True)
    counts = count_tp_fp(mt_tp_tmp, mt_fp_tmp)
    results['TP'] = counts[0]
    results['FP'] = counts[1]

    if var_type == 'snv':
        ratio = get_trans_untrans(mt_syn_tmp, pedigree, sample_list, mtdir)
        results['t_u_ratio'] = ratio

    return results


def filter_and_count(mt_tp: hl.MatrixTable, mt_fp: hl.MatrixTable, mt_syn: hl.MatrixTable, mt_roh: hl.MatrixTable, pedfile: str, mtdir: str) -> dict:
    '''
    Filter MT by various bins followed by genotype GQ and cauclate % of FP and TP remaining for each bifn
    :param hl.MatrixTable mt_tp: Input TP MatrixTable
    :param hl.MatrixTable mt_fp: Input FP MatrixTable
    :param hl.MatrixTable mt_syn: Input synonymous MatrixTable
    :param str mtdir: matrixtable directory
    :return: dict
    '''
    results = {'snv': {}, 'indel': {}}
    mtpath = mtdir.replace('file://', '')

    snp_mt_tp = mt_tp.filter_rows(mt_tp.type == snp_label)
    snp_mt_fp = mt_fp.filter_rows(mt_fp.type == snp_label)
    snp_mt_roh = mt_roh.filter_rows(mt_roh.type == snp_label)
    indel_mt_tp = mt_tp.filter_rows(mt_tp.type == indel_label)
    indel_mt_fp = mt_fp.filter_rows(mt_fp.type == indel_label)
    indel_mt_roh = mt_roh.filter_rows(mt_roh.type == indel_label)

    results['snv_total_tp'] = snp_mt_tp.count_rows()
    results['snv_total_fp'] = snp_mt_fp.count_rows()
    results['indel_total_tp'] = indel_mt_tp.count_rows()
    results['indel_total_fp'] = indel_mt_fp.count_rows()

    for bin in snp_bins:
        bin_str = "bin_" + str(bin)
        print("bin " + str(bin))
        mt_tp_tmp = snp_mt_tp.filter_rows(snp_mt_tp.info.rf_bin <= bin)
        tmpmtb1 = os.path.join(mtdir, "tmp1bx.mt")
        mt_tp_tmp = mt_tp_tmp.checkpoint(tmpmtb1, overwrite = True)

        mt_fp_tmp = snp_mt_fp.filter_rows(snp_mt_fp.info.rf_bin <= bin)
        tmpmtb2 = os.path.join(mtdir, "tmp2bx.mt")
        mt_fp_tmp = mt_fp_tmp.checkpoint(tmpmtb2, overwrite = True)

        mt_syn_tmp = mt_syn.filter_rows(mt_syn.info.rf_bin <= bin)
        tmpmtb3 = os.path.join(mtdir, "tmp3bx.mt")
        mt_syn_tmp = mt_syn_tmp.checkpoint(tmpmtb3, overwrite = True)

        mt_roh_tmp = snp_mt_roh.filter_rows(snp_mt_roh.info.rf_bin <= bin)
        tmpmtb4 = os.path.join(mtdir, "tmp4bx.mt")
        mt_roh_tmp = mt_roh_tmp.checkpoint(tmpmtb4, overwrite=True)

        for dp in dp_vals:
            dp_str = 'DP_' + str(dp)
            for gq in gq_vals:
                gq_str = 'GQ_' + str(gq)
                for ab in ab_vals:
                    ab_str = 'AB_' + str(ab)
                    print(dp_str + " " + gq_str + " " + ab_str)
                    filter_name = ("_").join([bin_str, dp_str, gq_str, ab_str])
                    snp_counts = filter_mt_count_tp_fp_t_u(mt_tp_tmp, mt_fp_tmp, mt_syn_tmp, pedfile, dp, gq, ab, 'snv', mtdir)
                    results['snv'][filter_name] = snp_counts

                    mt_roh_filtered = apply_hard_filters(mt_roh_tmp, dp=dp, gq=gq, ab=ab)
                    het_counts = count_hets_in_rohs(mt_roh_filtered, roh_path=roh_path)
                    write_out(het_counts, n_cores=240, out_prefix=os.path.join(mtdir, f'{filter_name}.snp.roh_stat'))

                    with open(os.path.join(mtpath, 'results.json'), 'w') as f:
                        json.dump(results, f)

    for bin in indel_bins:
        bin_str = "bin_" + str(bin)
        print("bin " + str(bin))

        mt_tp_tmp = indel_mt_tp.filter_rows(indel_mt_tp.info.rf_bin <= bin)
        tmpmtb1 = mtdir + "tmp1bx.mt"
        mt_tp_tmp = mt_tp_tmp.checkpoint(tmpmtb1, overwrite = True)

        mt_fp_tmp = indel_mt_fp.filter_rows(indel_mt_fp.info.rf_bin <= bin)
        tmpmtb2 = mtdir + "tmp2bx.mt"
        mt_fp_tmp = mt_fp_tmp.checkpoint(tmpmtb2, overwrite = True)

        mt_roh_tmp = indel_mt_roh.filter_rows(indel_mt_roh.info.rf_bin <= bin)
        tmpmtb4 = os.path.join(mtdir, "tmp4bx.mt")
        mt_roh_tmp = mt_roh_tmp.checkpoint(tmpmtb4, overwrite=True)

        for dp in dp_vals:
            dp_str = 'DP_' + str(dp)
            for gq in gq_vals:
                gq_str = 'GQ_' + str(gq)
                for ab in ab_vals:
                    ab_str = 'AB_' + str(ab)
                    print(dp_str + " " + gq_str + " " + ab_str)
                    filter_name = ("_").join([bin_str, dp_str, gq_str, ab_str])
                    indel_counts = filter_mt_count_tp_fp_t_u(mt_tp_tmp, mt_fp_tmp, mt_syn, pedfile, dp, gq, ab, 'indel', mtdir)
                    results['indel'][filter_name] = indel_counts

                    mt_roh_filtered = apply_hard_filters(mt_roh_tmp, dp=dp, gq=gq, ab=ab)
                    het_counts = count_hets_in_rohs(mt_roh_filtered, roh_path=roh_path)
                    write_out(het_counts, n_cores=240, out_prefix=os.path.join(mtdir, f'{filter_name}.indel.roh_stat'))

                    with open(os.path.join(mtpath, 'results.json'), 'w') as f:
                        json.dump(results, f)

    return results


def print_results(results: dict, outfile: str, vartype: str):
    '''
    Print results dict to a file
    :param dict results: results dict
    :param str outfile: output file path
    :param str vartype: variant type (snv or indel)
    '''
    header = ['filter', 'TP', 'FP']
    if vartype == 'snv':
        header = header + (['t_u_ratio'])

    with open(outfile, 'w') as o:
        o.write(("\t").join(header))
        o.write("\n")

        for var_f in results[vartype].keys():
            if vartype == 'snv':
                tp = str((results[vartype][var_f]['TP'] / results['snv_total_tp']) * 100)
                fp = str((results[vartype][var_f]['FP'] / results['snv_total_fp']) * 100)
            elif vartype == 'indel':
                tp = str((results[vartype][var_f]['TP'] / results['indel_total_tp']) * 100)
                fp = str((results[vartype][var_f]['FP'] / results['indel_total_fp']) * 100)
            outline = [var_f, tp, fp]
            if vartype == 'snv':
                tu = str(results[vartype][var_f]['t_u_ratio'])
                outline = outline + [tu]

            o.write(("\t").join(outline))
            o.write("\n")


def main():
    # set up
    inputs = parse_config()
    rf_dir = inputs['var_qc_rf_dir']
    mtdir = inputs['matrixtables_lustre_dir']
    resourcedir = inputs['resource_dir']
    annodir = inputs['annotation_lustre_dir']

    # initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    rf_htfile = os.path.join(rf_dir, rf_hash, "_gnomad_score_binning_tmp.ht")
    mtfile = os.path.join(mtdir, "mt_varqc_splitmulti.mt")
    cqfile = os.path.join(resourcedir, "all_consequences.txt")
    pedfile = 'file:///lustre/scratch123/projects/gnh_industry/Genes_and_Health_2023_02_44k/GH_44k_668-trios_QCed.mercury.fam'

    mt = hl.read_matrix_table(mtfile)
    mt = clean_mt(mt)
    mt_annot = annotate_with_rf(mt, rf_htfile)
    mt_annot = annotate_cq(mt_annot, cqfile)

    # mt_tp, mt_fp, mt_syn, mt_roh = filter_mts(mt_annot, roh_path, mtdir=os.path.join(mtdir, rf_hash))
    results = filter_and_count(
        # hl.read_matrix_table(mt_tp),
        # hl.read_matrix_table(mt_fp),
        # hl.read_matrix_table(mt_syn),
        # hl.read_matrix_table(mt_roh),
        hl.read_matrix_table(os.path.join(mtdir, rf_hash, 'tp.mt')),
        hl.read_matrix_table(os.path.join(mtdir, rf_hash, 'fp.mt')),
        hl.read_matrix_table(os.path.join(mtdir, rf_hash, 'syn.mt')),
        hl.read_matrix_table(os.path.join(mtdir, rf_hash, 'roh.mt')),
        pedfile, mtdir=os.path.join(mtdir, rf_hash)
    )

    outfile_snv = os.path.join(annodir, "genotype_hard_filter_comparison_snv.txt")
    outfile_indel = os.path.join(annodir, "genotype_hard_filter_comparison_indel.txt")
    print_results(results, outfile_snv, 'snv')
    print_results(results, outfile_indel, 'indel')


if __name__ == '__main__':
    main()
