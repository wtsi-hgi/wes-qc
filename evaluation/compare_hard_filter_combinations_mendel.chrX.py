import json
import os.path

import hail as hl
import pyspark
import datetime
from utils.utils import parse_config
from typing import Tuple, List, Dict, Any

from evaluation.compare_hard_filter_combinations import annotate_with_rf
from evaluation.compare_hard_filter_combinations import annotate_cq
from evaluation.compare_hard_filter_combinations import count_tp_fp
from utils.utils import rm_mt
from utils.utils import select_founders, collect_pedigree_samples
from utils.constants import Sex as sex

rf_hash = 'fa8798b1'
snp_bins = [80, 82, 83, 84]
indel_bins = [44, 46, 48, 49, 50]
gq_vals = [10, 15, 20]
dp_vals = [5, 10]
ab_vals = [0.2, 0.3]
missing_vals = [0, 0.5, 0.9, 0.95]
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


def filter_mts(mt: hl.MatrixTable, mtdir: str) -> Tuple[hl.MatrixTable, hl.MatrixTable]:
    """
    Split matrixtable and return tables with just TP, just FP, just synonymous
    :param hl.MatrixTable mt: Input mtfile
    :param st mtdir: matrixtable directory
    :return: tuple of 4 hl.MatrixTable objects
    """
    cores_in_cluster = 236

    now = datetime.datetime.now()
    print(now.time())

    mt_true = mt.filter_rows(mt.TP == True)  # TP variants
    mt_true = hl.variant_qc(mt_true)  # remove unused rows
    mt_true = mt_true.filter_rows(mt_true.variant_qc.n_non_ref == 0, keep=False)
    mt_true = mt_true.drop(mt_true.variant_qc)

    mt_false = mt.filter_rows(mt.FP == True)  # FP variants
    mt_false = hl.variant_qc(mt_false)  # remove unused rows
    mt_false = mt_false.filter_rows(mt_false.variant_qc.n_non_ref == 0, keep=False)
    mt_false = mt_false.drop(mt_false.variant_qc)

    tmpmtt = os.path.join(mtdir, "tp.mt")
    tmpmtf = os.path.join(mtdir, "fp.mt")

    mt_true = mt_true.naive_coalesce(cores_in_cluster).checkpoint(tmpmtt, overwrite=True)
    mt_false = mt_false.naive_coalesce(cores_in_cluster).checkpoint(tmpmtf, overwrite=True)

    return mt_true, mt_false


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


def get_count_tp_fp(mt_tp: hl.MatrixTable, mt_fp: hl.MatrixTable):
    '''
    Filter mt by each rf bin in a list, then genotype hard filters and count remaining TP and P variants
    :param hl.MatrixTable mt_tp: Input TP MatrixTable
    :param hl.MatrixTable mt_fp: Input FP MatrixTable
    :return: Dict containing bin and remaining TP/FP count
    '''
    results = {}

    counts = count_tp_fp(mt_tp, mt_fp)
    results['TP'] = counts[0]
    results['FP'] = counts[1]

    return results


def get_trans_untrans(mt: hl.MatrixTable, pedigree: hl.Pedigree, mtdir: str, gender: str) -> float:
    """
    get transmitted/untransmitted ratio
    :param hl.MatrixTable mt: matrixtable
    :param hl.Pedigree pedigree: Hail Pedigree
    :param str mtdir: matrixtable directory
    :return float:
    """
    parts = min(mt.n_partitions(), 2000)

    # filter to synonymous
    mt_syn = mt.filter_rows(mt.consequence == 'synonymous_variant')  # synonymous for transmitted/unstransmitted
    mt_syn = hl.variant_qc(mt_syn)  # remove unused rows
    mt_syn = mt_syn.filter_rows(mt_syn.variant_qc.n_non_ref == 0, keep=False)

    # restrict to samples in trios, annotate with AC and filter to AC == 1 in parents
    sample_list = collect_pedigree_samples(pedigree)  # list of samples in trios
    mt2 = mt_syn.filter_cols(hl.set(sample_list).contains(mt_syn.s))

    founders = select_founders(pedigree)
    mt_founders = mt2.filter_cols(hl.set(founders).contains(mt2.s))
    mt_founders = hl.variant_qc(mt_founders, name='varqc_founders')

    mt2 = mt2.annotate_rows(varqc_trios=hl.Struct(AC=mt_founders.index_rows(mt2.row_key).varqc_founders.AC))

    # split to potentially transmitted/untransmitted
    trans_mt = mt2.filter_rows(mt2.varqc_trios.AC[1] == 1)
    tmpmt5 = os.path.join(mtdir, "tmp5x.mt")
    trans_mt = trans_mt.checkpoint(tmpmt5, overwrite=True)

    # run tdt function for potential trans and untrans
    tdt_ht = hl.transmission_disequilibrium_test(trans_mt, pedigree)

    gender_mt = trans_mt.filter_cols(trans_mt.sex == gender)
    melted_mt = gender_mt.entries().select('GT')

    # wrong
    # trans_ht = melted_mt.filter(melted_mt.GT.is_hom_ref(), keep=False).key_by('locus', 'alleles').join(tdt_ht).distinct()
    # untrans_ht = melted_mt.filter(melted_mt.GT.is_hom_ref(), keep=True).key_by('locus', 'alleles').join(tdt_ht).distinct()

    trans = trans_ht.aggregate(hl.agg.sum(trans_ht.t))
    untrans = untrans_ht.aggregate(hl.agg.sum(untrans_ht.u))
    if untrans > 0:
        ratio = trans/untrans
    else:
        ratio = 0

    return ratio


def walk_grid(mt: hl.MatrixTable, bins: List[int], results: Dict[str, Any], pedigree: hl.Pedigree,
              wd: str, outjson: str, mutation: str, resume=False):
    assert mutation in ('snv', 'indel')
    mt_snp_hard_path = os.path.join(wd, 'tmp.hard_filters_combs.snp-hard.mt')
    genders = (sex.male, sex.female)

    for bin in bins:
        bin_str = "bin_" + str(bin)
        print("bin " + str(bin))
        mt_snp_bin = mt.filter_rows(mt.info.rf_bin <= bin)

        for dp in dp_vals:
            dp_str = f'DP_{dp}'
            for gq in gq_vals:
                gq_str = f'GQ_{gq}'
                for ab in ab_vals:
                    ab_str = f'AB_{ab}'
                    for missing in missing_vals:
                        missing_str = f'missing_{missing}'
                        print(f'{dp_str} {gq_str} {ab_str} {missing_str}')

                        filter_prefix = "_".join([bin_str, dp_str, gq_str, ab_str, missing_str])

                        if resume:
                            filter_results = [x for x in results[mutation].keys() if x.startswith(filter_prefix)]
                            if len(filter_results) == len(genders):
                                continue

                        mt_snp_hard = apply_hard_filters(mt_snp_bin, dp=dp, gq=gq, ab=ab, call_rate=missing)
                        mt_snp_hard = mt_snp_hard.checkpoint(mt_snp_hard_path, overwrite=True)

                        _, _, per_sample, _ = hl.mendel_errors(mt_snp_hard.GT, pedigree)
                        per_sample = per_sample.join(mt_snp_hard.cols().select('sex'))

                        for gender in genders:
                            sex_str = f'sex_{gender}'
                            filter_name = f'{filter_prefix}_{sex_str}'

                            # if mutation == 'snv':
                            #     ratio = get_trans_untrans(mt_snp_hard, pedigree, wd, gender=gender)
                            #     results[mutation][filter_name]['t_u_ratio'] = ratio

                            mt_snp_sex = mt_snp_hard.filter_cols(mt_snp_hard.sex == gender)
                            mt_tp_tmp, mt_fp_tmp = filter_mts(mt_snp_sex, mtdir=wd)

                            results[mutation][filter_name] = get_count_tp_fp(mt_tp_tmp, mt_fp_tmp)

                            if gender == sex.male:
                                results[mutation][filter_name]['hets_fraction'] = count_male_hets(mt_snp_sex)

                            per_sample_gender = per_sample.filter(per_sample.sex == gender)
                            mendel_errors = per_sample_gender.aggregate(hl.agg.sum(per_sample_gender.errors))
                            results[mutation][filter_name]['mendel_errors'] = mendel_errors

                            with open(outjson, 'w') as f:
                                json.dump(results, f)

    rm_mt(mt_snp_hard_path)
    return results


def count_male_hets(mt: hl.MatrixTable) -> float:
    mt_male = mt.filter_cols(mt.sex == sex.male)
    mt_nonpar = mt_male.filter_rows(mt_male.locus.in_x_nonpar())
    mt_nonpar = mt_nonpar.annotate_cols(hets_fraction=hl.agg.fraction(mt_nonpar.GT.is_het()))
    mean_hets_fraction = mt_nonpar.aggregate_cols(hl.agg.mean(mt_nonpar.hets_fraction))
    return mean_hets_fraction


def filter_and_count(mt_path: str, pedfile: str, mtdir: str, resume=False) -> dict:
    """
    Filter MT by various bins followed by genotype GQ and calculate % of FP and TP remaining for each bin
    :return: dict
    """
    mtpath = mtdir.replace('file://', '')
    out_json = os.path.join(mtpath, 'results-chrX.json')

    if resume:
        with open(out_json, 'r') as f:
            results = json.load(f)
    else:
        results = {'snv': {}, 'indel': {}}

    mt = hl.read_matrix_table(mt_path)
    pedigree = hl.Pedigree.read(pedfile)

    mt_snp = mt.filter_rows(mt.type == snp_label)
    mt_indel = mt.filter_rows(mt.type == indel_label)

    snp_mt_tp, snp_mt_fp = filter_mts(mt_snp, mtdir=mtdir)
    results['snv_total_tp'] = snp_mt_tp.count_rows()
    results['snv_total_fp'] = snp_mt_fp.count_rows()

    results = walk_grid(
        mt=mt_snp, bins=snp_bins, mutation='snv',
        pedigree=pedigree, results=results,
        wd=mtdir, outjson=out_json, resume=resume
    )

    indel_mt_tp, indel_mt_fp = filter_mts(mt_indel, mtdir=mtdir)
    results['indel_total_tp'] = indel_mt_tp.count_rows()
    results['indel_total_fp'] = indel_mt_fp.count_rows()

    results = walk_grid(
        mt=mt_indel,  bins=indel_bins, mutation='indel',
        pedigree=pedigree, results=results,
        wd=mtdir, outjson=out_json
    )

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


def read_genders(sexfile: str, linkfile: str) -> hl.Table:
    genders = hl.import_fam(sexfile).key_by('fam_id')
    link = hl.import_table(linkfile, key='OrageneID.GSA')
    genders = genders.join(link).key_by('OrageneID.WES')
    return genders


def annotate_sex(mt: hl.MatrixTable, genders: hl.Table) -> hl.MatrixTable:
    sex_expr = (hl.switch(genders[mt.s].is_female)
                  .when(True, sex.female)
                  .when(False, sex.male)
                  .when_missing('unknown')
                  .default('error'))
    mt = mt.annotate_cols(sex=sex_expr)
    return mt


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
    mtfile = os.path.join(mtdir, "mt_varqc_splitmulti.chrX.mt")
    cqfile = os.path.join(resourcedir, "all_consequences.txt")
    pedfile = 'file:///lustre/scratch123/projects/gnh_industry/Genes_and_Health_2023_02_44k/GH_44k_668-trios_QCed.mercury.consistent.fam'
    wd = os.path.join(mtdir, rf_hash)

    sexfile = 'file:///lustre/scratch123/projects/gnh_industry/Genes_and_Health_2023_02_44k/samples_sex.fam'
    linkfile = 'file:///lustre/scratch123/projects/gnh_industry/Genes_and_Health_2023_02_44k/link_OrageneID_all-WES_GSA.txt'
    genders = read_genders(sexfile, linkfile)

    mt = hl.read_matrix_table(mtfile)
    mt = clean_mt(mt)
    mt_annot = annotate_with_rf(mt, rf_htfile)
    mt_annot = annotate_cq(mt_annot, cqfile)
    mt_annot = annotate_sex(mt_annot, genders)

    mt_annot_path = os.path.join(wd, 'tmp.hard_filters_combs.mt')
    # mt_annot.write(mt_annot_path, overwrite=True)

    results = filter_and_count(
        mt_path=mt_annot_path,
        pedfile=pedfile,
        mtdir=wd,
        resume=True
    )

    outfile_snv = os.path.join(annodir, "genotype_hard_filter_comparison_snv.txt")
    outfile_indel = os.path.join(annodir, "genotype_hard_filter_comparison_indel.txt")
    # print_results(results, outfile_snv, 'snv')
    # print_results(results, outfile_indel, 'indel')


if __name__ == '__main__':
    main()
