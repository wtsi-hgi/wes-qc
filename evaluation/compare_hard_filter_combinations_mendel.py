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
from utils.utils import rm_mt

rf_hash = 'fa8798b1'
snp_bins = [80, 82, 83, 84]
indel_bins = [44, 46, 48, 49, 50]
gq_vals = [10, 15, 20]
dp_vals = [5, 10]
ab_vals = [0.2, 0.3]
missing_vals = [0.9, 0.95]
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


def filter_mts(mt: hl.MatrixTable, roh: hl.Table, mtdir: str) -> Tuple[hl.MatrixTable, hl.MatrixTable, hl.MatrixTable, hl.MatrixTable]:
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

    mt_true = mt_true.checkpoint(tmpmtt, overwrite=True)
    mt_false = mt_false.checkpoint(tmpmtf, overwrite=True)
    mt_syn = mt_syn.checkpoint(tmpmts, overwrite=True)

    samples = roh.aggregate(hl.agg.collect_as_set(roh.ID))
    mt_roh = mt.filter_cols(hl.set(samples).contains(mt.s))
    mt_roh = mt_roh.checkpoint(tmpmtr, overwrite=True)

    return mt_true, mt_false, mt_syn, mt_roh


def count_rows(mt: hl.MatrixTable, mutation_type: str) -> int:
    c = mt.filter_rows(mt.type == mutation_type).count_rows()
    return c


def apply_hard_filters(mt: hl.MatrixTable, dp: int, gq: int, ab: float, call_rate: float) -> hl.MatrixTable:
    filter_condition = (
            (mt.GT.is_het() & (mt.HetAB < ab)) |
            (mt.DP < dp) |
            (mt.GQ < gq)
    )
    mt_tmp = mt.annotate_entries(
        hard_filters=hl.if_else(filter_condition, 'Fail', 'Pass')
    )
    mt_tmp = mt_tmp.filter_entries(mt_tmp.hard_filters == 'Pass')

    mt_tmp = mt_tmp.annotate_rows(pass_count=hl.agg.count_where(mt_tmp.hard_filters == 'Pass'))
    mt_tmp = mt_tmp.filter_rows(mt_tmp.pass_count/mt_tmp.count_cols() > call_rate)

    return mt_tmp


def filter_mt_count_tp_fp_t_u(mt_tp: hl.MatrixTable, mt_fp: hl.MatrixTable, mt_syn: hl.MatrixTable, pedigree: hl.Pedigree, var_type: str, mtdir: str):
    '''
    Filter mt by each rf bin in a list, then genotype hard filters and count remaining TP and P variants
    :param hl.MatrixTable mt_tp: Input TP MatrixTable
    :param hl.MatrixTable mt_fp: Input FP MatrixTable
    :param hl.MatrixTable mt_syn: Input synonymous MatrixTable
    :param str pedigree: hail pedigree object
    :param int dp: DP threshold
    ;param int gq; GQ threshold
    :param float ab: allele balance threshold
    :param  str var_type: variant type (snv/indel)
    :param str mtdir: matrixtable directory
    :return: Dict containing bin and remaning TP/FP count
    '''
    results = {}

    #genotype hard filters - should put this in a different method
    now = datetime.datetime.now()
    print(now.time())

        #remove unused rows
    mt_tp_tmp = hl.variant_qc(mt_tp)
    mt_tp_tmp = mt_tp_tmp.filter_rows(mt_tp_tmp.variant_qc.n_non_ref == 0, keep = False)

        #remove unused rows
    mt_fp_tmp = hl.variant_qc(mt_fp)
    mt_fp_tmp = mt_fp_tmp.filter_rows(mt_fp_tmp.variant_qc.n_non_ref == 0, keep = False)

        #remove unused rows
    mt_syn_tmp = hl.variant_qc(mt_syn)
    mt_syn_tmp = mt_syn_tmp.filter_rows(mt_syn_tmp.variant_qc.n_non_ref == 0, keep = False)

    counts = count_tp_fp(mt_tp_tmp, mt_fp_tmp)
    results['TP'] = counts[0]
    results['FP'] = counts[1]

    if var_type == 'snv':
        ratio = get_trans_untrans(mt_syn_tmp, pedigree, mtdir)
        results['t_u_ratio'] = ratio

    return results


def filter_and_count(mt_path: str, roh_path: str, pedfile: str, mtdir: str) -> dict:
    '''
    Filter MT by various bins followed by genotype GQ and cauclate % of FP and TP remaining for each bin
    :param hl.MatrixTable mt_tp: Input TP MatrixTable
    :param hl.MatrixTable mt_fp: Input FP MatrixTable
    :param hl.MatrixTable mt_syn: Input synonymous MatrixTable
    :param str mtdir: matrixtable directory
    :return: dict
    '''
    results = {'snv': {}, 'indel': {}}
    mtpath = mtdir.replace('file://', '')

    mt = hl.read_matrix_table(mt_path)
    roh = import_roh(roh_path)
    pedigree = hl.Pedigree.read(pedfile)

    mt_snp = mt.filter_rows(mt.type == snp_label)
    mt_indel = mt.filter_rows(mt.type == indel_label)

    mt_snp_path = os.path.join(mtdir, 'tmp.hard_filters_combs.snp.mt')
    mt_snp = mt_snp.checkpoint(mt_snp_path, overwrite=True)
    mt_indel_path = os.path.join(mtdir, 'tmp.hard_filters_combs.indel.mt')
    mt_indel = mt_indel.repartition(mt.n_partitions() // 4).checkpoint(mt_indel_path, overwrite=True)

    snp_mt_tp, snp_mt_fp, _, _ = filter_mts(mt_snp, roh, mtdir=mtdir)
    results['snv_total_tp'] = snp_mt_tp.count_rows()
    results['snv_total_fp'] = snp_mt_fp.count_rows()

    for bin in snp_bins:
        bin_str = "bin_" + str(bin)
        print("bin " + str(bin))
        mt_snp_bin = mt_snp.filter_rows(mt_snp.info.rf_bin <= bin)

        for dp in dp_vals:
            dp_str = 'DP_' + str(dp)
            for gq in gq_vals:
                gq_str = 'GQ_' + str(gq)
                for ab in ab_vals:
                    ab_str = 'AB_' + str(ab)
                    for missing in missing_vals:
                        missing_str = f'missing_{missing}'
                        print(f'{dp_str} {gq_str} {ab_str} {missing_str}')
                        filter_name = ("_").join([bin_str, dp_str, gq_str, ab_str, missing_str])

                        mt_snp_hard = apply_hard_filters(mt_snp_bin, dp=dp, gq=gq, ab=ab, call_rate=missing)
                        mt_snp_hard_path = os.path.join(mtdir, 'tmp.hard_filters_combs.snp-hard.mt')
                        mt_snp_hard = mt_snp_hard.checkpoint(mt_snp_hard_path, overwrite=True)

                        mt_tp_tmp, mt_fp_tmp, mt_syn_tmp, mt_roh_tmp = filter_mts(mt_snp_hard, roh, mtdir=mtdir)

                        snp_counts = filter_mt_count_tp_fp_t_u(mt_tp_tmp, mt_fp_tmp, mt_syn_tmp, pedigree, 'snv', mtdir)
                        results['snv'][filter_name] = snp_counts

                        # if not os.path.exists(os.path.join(mtpath, f'{filter_name}.snp.roh_stat.tsv')):
                        #     mt_roh_filtered = apply_hard_filters(mt_roh_tmp, dp=dp, gq=gq, ab=ab)
                        #     het_counts = count_hets_in_rohs(mt_roh_filtered, roh_path=roh_path)
                        #     write_out(het_counts, n_cores=240, out_prefix=os.path.join(mtdir, f'{filter_name}.snp.roh_stat'))

                        all_errors, _, _, _ = hl.mendel_errors(mt_snp_hard.GT, pedigree)
                        results['snv'][filter_name]['mendel_errors'] = all_errors.count()

                        with open(os.path.join(mtpath, 'results-mendel-missing.json'), 'w') as f:
                            json.dump(results, f)

    indel_mt_tp, indel_mt_fp, _, _ = filter_mts(mt_indel, roh, mtdir=mtdir)
    results['indel_total_tp'] = indel_mt_tp.count_rows()
    results['indel_total_fp'] = indel_mt_fp.count_rows()

    for bin in indel_bins:
        bin_str = "bin_" + str(bin)
        print("bin " + str(bin))
        mt_indel_bin = mt_indel.filter_rows(mt_indel.info.rf_bin <= bin)

        for dp in dp_vals:
            dp_str = 'DP_' + str(dp)
            for gq in gq_vals:
                gq_str = 'GQ_' + str(gq)
                for ab in ab_vals:
                    ab_str = 'AB_' + str(ab)
                    for missing in missing_vals:
                        missing_str = f'missing_{missing}'
                        print(f'{dp_str} {gq_str} {ab_str} {missing_str}')
                        filter_name = ("_").join([bin_str, dp_str, gq_str, ab_str, missing_str])

                        mt_indel_hard = apply_hard_filters(mt_indel_bin, dp=dp, gq=gq, ab=ab, call_rate=missing)
                        mt_indel_hard_path = os.path.join(mtdir, 'tmp.hard_filters_combs.indel-hard.mt')
                        mt_indel_hard = mt_indel_hard.checkpoint(mt_indel_hard_path, overwrite=True)

                        mt_tp_tmp, mt_fp_tmp, mt_syn_tmp, mt_roh_tmp = filter_mts(mt_indel_hard, roh, mtdir=mtdir)

                        indel_counts = filter_mt_count_tp_fp_t_u(mt_tp_tmp, mt_fp_tmp, mt_syn_tmp, pedigree, 'indel', mtdir)
                        results['indel'][filter_name] = indel_counts

                        # if not os.path.exists(os.path.join(mtpath, f'{filter_name}.indel.roh_stat.tsv')):
                        #     mt_roh_filtered = apply_hard_filters(mt_roh_tmp, dp=dp, gq=gq, ab=ab)
                        #     het_counts = count_hets_in_rohs(mt_roh_filtered, roh_path=roh_path)
                        #     write_out(het_counts, n_cores=240, out_prefix=os.path.join(mtdir, f'{filter_name}.indel.roh_stat'))

                        all_errors, _, _, _ = hl.mendel_errors(mt_indel_hard.GT, pedigree)
                        results['indel'][filter_name]['mendel_errors'] = all_errors.count()

                        with open(os.path.join(mtpath, 'results-mendel-missing.json'), 'w') as f:
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
    annodir = inputs['annotation_lustre_dir_local']

    # initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    rf_htfile = os.path.join(rf_dir, rf_hash, "_gnomad_score_binning_tmp.ht")
    mtfile = os.path.join(mtdir, "mt_varqc_splitmulti.mt")
    cqfile = os.path.join(resourcedir, "all_consequences.txt")
    pedfile = 'file:///lustre/scratch123/projects/gnh_industry/Genes_and_Health_2023_02_44k/GH_44k_668-trios_QCed.mercury.consistent.fam'
    wd = os.path.join(mtdir, rf_hash)

    mt = hl.read_matrix_table(mtfile)
    mt = clean_mt(mt)
    mt_annot = annotate_with_rf(mt, rf_htfile)
    mt_annot = annotate_cq(mt_annot, cqfile)

    mt_annot_path = os.path.join(wd, 'tmp.hard_filters_combs.mt')
    # mt_annot.write(mt_annot_path, overwrite=True)

    results = filter_and_count(
        mt_path=mt_annot_path,
        roh_path=roh_path,
        pedfile=pedfile,
        mtdir=wd
    )

    outfile_snv = os.path.join(annodir, "genotype_hard_filter_comparison_snv.txt")
    outfile_indel = os.path.join(annodir, "genotype_hard_filter_comparison_indel.txt")
    # print_results(results, outfile_snv, 'snv')
    # print_results(results, outfile_indel, 'indel')


if __name__ == '__main__':
    main()
