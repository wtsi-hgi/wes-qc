# generate truth sets for variant QC random forest
import hail as hl
import pyspark
import argparse
from typing import Tuple
from utils.utils import parse_config, rm_mt
from gnomad.utils.filtering import filter_to_adj, filter_to_autosomes
from gnomad.utils.annotations import unphase_call_expr, bi_allelic_site_inbreeding_expr, add_variant_type, annotate_adj, bi_allelic_expr
from gnomad.sample_qc.relatedness import generate_trio_stats_expr


def get_options():
    '''
    Get options from the command line
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("--truth",
        help="create truth set HT", action="store_true")
    parser.add_argument("--annotation",
        help="family, inbreeding, allele data and dnm annotation", action="store_true")
    parser.add_argument("--all",
        help="run all steps", action="store_true")
    args = parser.parse_args()
    if not args.all and not args.truth and not args.annotation:
        print("No args specified, either --truth, --annotation or --all must be used")
        exit(1)

    return args


def get_truth_ht(omni, mills, thousand_genomes, hapmap, truth_ht_file):
    '''
    Generate truth set HT
    :param str omni: file containing kgp_omni 1000 Genomes intersection Onni 2.5M array
    :param str mills: file containing Mills & Devine indels
    :param str thousand_genomes: file containing high confidence sites in 1000 genonmes
    :param str hapmap: file containing hapmap
    :param str truth_ht_file: output file name
    '''
    omni_ht = hl.read_table(omni)
    mills_ht = hl.read_table(mills)
    thousand_genomes_ht = hl.read_table(thousand_genomes)
    hapmap_ht = hl.read_table(hapmap)

    truth_ht = hapmap_ht.select(hapmap=True).join(
        omni_ht.select(omni=True), how="outer").join(
            thousand_genomes_ht.select(kgp_phase1_hc=True), how="outer").join(
                mills_ht.select(mills=True), how="outer").repartition(200, shuffle=False)
    truth_ht.write(truth_ht_file, overwrite=True)
    

def split_multi_and_var_qc(mtfile: str, varqc_mtfile: str, varqc_mtfile_split: str):
    '''
    Adapted from https://github.com/broadinstitute/gnomad_qc/blob/3d79bdf0f7049c209b4659ff8c418a1b859d7cfa/gnomad_qc/v2/annotations/generate_qc_annotations.py
    Run split_multi_hts and variant QC on inout mt after sample QC
    :param str mtfile: Input mtfile, raw variants with sample QC fails removed
    :param str varqc_mtfile: Output mt with adj annotation 
    :param str varqc_mtfile_split: Output mt with variant QC annotation and split multiallelics
    '''
    mt = hl.read_matrix_table(mtfile)

    #        # restrict to samples in trios
    # samplefile = "file:///lustre/scratch123/qc/resources/samples_in_trios.txt"
    # sampleht = hl.import_table(samplefile)
    # sampleht = sampleht.key_by('s')
    # mt = mt.filter_cols(hl.is_defined(sampleht[mt.s]))
    
    #before splitting annotate with sum_ad
    mt = mt.annotate_entries(sum_AD = hl.sum(mt.AD))
    mt = annotate_adj(mt)
    mt = mt.checkpoint(varqc_mtfile, overwrite=True)
    mt = hl.split_multi_hts(mt)
    
    tmp_mt = varqc_mtfile_split + "_tmp"
    print("writing split mt")
    mt = mt.checkpoint(tmp_mt, overwrite=True)

    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(mt.variant_qc.n_non_ref == 0, keep = False) 

    # mt = hl.split_multi_hts(mt)
    # tmp_mt = varqc_mtfile_split + "_tmp"
    # print("writing split mt")
    # mt.write(tmp_mt, overwrite=True)

    # #add AC for trtios only
    # mt = mt.annotate_rows(trioAC = hl.agg.sum(mt.AD[1]))

    #min(GT_AD) used here - but could change this to cover just those with few alts - moved to after split
    print("Annotating entries with allele balance")
    mt = mt.annotate_entries(
        HetAB=hl.case()
                .when(mt.sum_AD == 0, hl.missing('float64'))
                .when(mt.GT.is_het(), hl.min(mt.AD[1]) / mt.sum_AD)  # hl.min(mt.AD[1]) / hl.sum(mt.AD)
                .or_missing()
    )

    print("Annotating variants with mean allele balance")
    mt = mt.annotate_rows(
        meanHetAB=hl.agg.mean(mt.HetAB)
    )

    # replacing NaN values with NA
    mt = mt.annotate_rows(
        meanHetAB=hl.if_else(hl.is_nan(mt.meanHetAB), hl.missing('float64'), mt.meanHetAB)
    )

    print("writing split mt")
    mt.write(varqc_mtfile_split, overwrite=True)
    rm_mt(tmp_mt)


def read_fam(fam_file: str) -> hl.Table:
    '''
    Taken from https://github.com/broadinstitute/gnomad_qc/blob/3d79bdf0f7049c209b4659ff8c418a1b859d7cfa/gnomad_qc/v2/annotations/generate_qc_annotations.py
    '''
    columns = ['fam_id', 's', 'pat_id', 'mat_id', 'is_female']
    return hl.import_table(fam_file, no_header=True).rename({f'f{i}': c for i, c in enumerate(columns)}).key_by('s')


def annotate_unrelated_sample(mt: hl.MatrixTable, fam_file: str) -> hl.MatrixTable:
    '''
    Taken from https://github.com/broadinstitute/gnomad_qc/blob/3d79bdf0f7049c209b4659ff8c418a1b859d7cfa/gnomad_qc/v2/annotations/generate_qc_annotations.py
    '''
    fam_ht = read_fam(fam_file)
    return mt.annotate_cols(unrelated_sample=hl.is_missing(fam_ht[mt.s]))


def family_stats(mt: hl.MatrixTable, ped: hl.Pedigree, group_name: str) -> Tuple[hl.expr.StructExpression, hl.Table]:
    '''
    Taken from https://github.com/broadinstitute/gnomad_qc/blob/3d79bdf0f7049c209b4659ff8c418a1b859d7cfa/gnomad_qc/v2/annotations/generate_qc_annotations.py
    '''
    tdt_table = hl.transmission_disequilibrium_test(mt, ped)
    _, _, per_sample, per_variant = hl.mendel_errors(mt.GT, ped)
    family_stats_struct = hl.struct(mendel=per_variant[mt.row_key],
                                    tdt=tdt_table[mt.row_key],
                                    unrelated_qc_callstats=hl.agg.filter(mt.unrelated_sample,
                                                                         hl.agg.call_stats(mt.GT, mt.alleles)),
                                    meta={'group': group_name})
    return family_stats_struct, per_sample


def generate_family_stats(mt: hl.MatrixTable, fam_file: str, calculate_adj: bool = False) -> Tuple[hl.Table, hl.Table]:
    """
    Adapted from https://github.com/broadinstitute/gnomad_qc/blob/3d79bdf0f7049c209b4659ff8c418a1b859d7cfa/gnomad_qc/v2/annotations/generate_qc_annotations.py
    Writes bi-allelic sites MT with the following annotations:
     - family_stats (TDT, Mendel Errors, AC_unrelated_qc)
     - truth_data (presence in Omni, HapMap, 1KG/TGP high conf SNVs, Mills)
    :param MatrixTable mt: Full MT
    :param str fam_file: Fam pedigree file location
    :param bool calculate_adj: Whether to also calculate family metrics for adj genotypes
    :return: Table with qc annotations
    :rtype: Table
    """
    #mt = mt.select_cols(high_quality=mt.meta.high_quality)
    mt = mt.select_rows()
    mt = annotate_unrelated_sample(mt, fam_file)

    # Unphased for now, since mendel_errors does not support phased alleles
    mt = mt.annotate_entries(GT=unphase_call_expr(mt.GT))
    ped = hl.Pedigree.read(fam_file, delimiter='\\t')
    family_stats_struct, family_stats_sample_ht = family_stats(mt, ped, 'raw')
    mt = mt.annotate_rows(family_stats=[family_stats_struct])

    if calculate_adj:
        mt = filter_to_adj(mt)
        adj_family_stats_struct, adj_family_stats_sample_ht = family_stats(mt, ped, 'adj')

        family_stats_sample_ht = family_stats_sample_ht.annotate(
            adj=adj_family_stats_sample_ht[family_stats_sample_ht.s])

        mt = mt.annotate_rows(family_stats=mt.family_stats.append(adj_family_stats_struct))

    return mt.rows(), family_stats_sample_ht


def generate_trio_stats(
    mt: hl.MatrixTable, autosomes_only: bool = True, bi_allelic_only: bool = True
) -> hl.Table:
    """
    Adapted from https://github.com/broadinstitute/gnomad_qc/blob/3d79bdf0f7049c209b4659ff8c418a1b859d7cfa/gnomad_qc/v2/annotations/generate_qc_annotations.py
    Default function to run `generate_trio_stats_expr` to get trio stats stratified by raw and adj
    .. note::
        Expects that `mt` is it a trio matrix table that was annotated with adj and if dealing with
        a sparse MT `hl.experimental.densify` must be run first.
        By default this pipeline function will filter `mt` to only autosomes and bi-allelic sites.
    :param mt: A Trio Matrix Table returned from `hl.trio_matrix`. Must be dense
    :param autosomes_only: If set, only autosomal intervals are used.
    :param bi_allelic_only: If set, only bi-allelic sites are used for the computation
    :return: Table with trio stats
    """
    if autosomes_only:
        mt = mt.filter_rows(mt.locus.in_autosome())
    #    mt = filter_to_autosomes(mt)
    if bi_allelic_only:
        mt = mt.filter_rows(bi_allelic_expr(mt))

    trio_adj = mt.proband_entry.adj & mt.father_entry.adj & mt.mother_entry.adj

    ht = mt.select_rows(
        **generate_trio_stats_expr(
            mt,
            transmitted_strata={"raw": True, "adj": trio_adj},
            de_novo_strata={"raw": True, "adj": trio_adj},
            ac_strata={"raw": True, "adj": trio_adj},
        )
    ).rows()

    return ht


def trio_family_dnm_annotation(varqc_mtfile: str, pedfile: str, trio_mtfile: str, trio_stats_htfile: str, fam_stats_htfile: str, fam_stats_mtfile: str, fam_stats_gnomad_mtfile: str, gnomad_htfile: str, dnm_htfile: str):
    '''
    Create matrixtable of just trios, create family stats
    :param str varqc_mtfile: Input matrixtable file
    :param str pedfile: Plink-format ped file
    :param str trio_mtfile: Trio matrixtable file
    :param str trio_stats_htfile: Trio stats hail table file
    :param str fam_stats_htfile: Family stats hail table file
    :param str fam_stats_mtfile: Family stats matrixtble file
    :param str fam_stats_gnomad_mtfile: Family statts with gnomad annotation matrixtable file
    :param str gnomad_htfile: Gnomad AF hail table
    :param str dnm_htfile: De novo hail table file
    '''
    pedigree = hl.Pedigree.read(pedfile)
    mt = hl.read_matrix_table(varqc_mtfile)

    #filter mt to samples that are in trios only and re-run varqc
    trio_sample_ht = hl.import_fam(pedfile)
    sample_list = trio_sample_ht.id.collect() + trio_sample_ht.pat_id.collect() + trio_sample_ht.mat_id.collect()
    mt2 = mt.filter_cols(hl.set(sample_list).contains(mt.s))
    mt2 = hl.variant_qc(mt2, name = 'varqc_trios')

    trio_dataset = hl.trio_matrix(mt2, pedigree, complete_trios=True)
    trio_dataset.write(trio_mtfile, overwrite=True)

    trio_stats_ht = generate_trio_stats(trio_dataset, autosomes_only=True, bi_allelic_only=False)
    trio_stats_ht.write(trio_stats_htfile, overwrite=True)

    print("Generating family stats")
    (ht1, famstats_ht) = generate_family_stats(mt, pedfile)
    ht1.write(fam_stats_htfile ,overwrite=True)

    mt = mt.annotate_rows(family_stats=ht1[mt.row_key].family_stats)
    mt = mt.checkpoint(fam_stats_mtfile, overwrite=True)
    #add gnomad AFs
    gnomad_ht = hl.read_table(gnomad_htfile)
    mt = mt.annotate_rows(gnomad_maf=gnomad_ht[mt.row_key].freq[0].AF)
    mt.write(fam_stats_gnomad_mtfile, overwrite=True)
    #make DNM table
    de_novo_table = hl.de_novo(mt, pedigree, mt.gnomad_maf)
    de_novo_table = de_novo_table.key_by('locus', 'alleles').collect_by_key('de_novo_data')
    de_novo_table.repartition(480).write(dnm_htfile, overwrite=True)

    rm_mt(fam_stats_mtfile)


def generate_allele_data(mt: hl.MatrixTable) -> hl.Table:
    """
    Taken from https://github.com/broadinstitute/gnomad_qc/blob/3d79bdf0f7049c209b4659ff8c418a1b859d7cfa/gnomad_qc/v2/annotations/generate_qc_annotations.py
    Writes bi-allelic sites MT with the following annotations:
     - allele_data (nonsplit_alleles, has_star, variant_type, and n_alt_alleles)
    :param MatrixTable mt: Full unsplit MT
    :return: Table with allele data annotations
    :rtype: Table
    """
    ht = mt.rows().select()
    allele_data = hl.struct(nonsplit_alleles=ht.alleles,
                            has_star=hl.any(lambda a: a == '*', ht.alleles))
    ht = ht.annotate(allele_data=allele_data.annotate(
        **add_variant_type(ht.alleles)))

    ht = hl.split_multi_hts(ht)
    allele_type = (hl.case()
                   .when(hl.is_snp(ht.alleles[0], ht.alleles[1]), 'snv')
                   .when(hl.is_insertion(ht.alleles[0], ht.alleles[1]), 'ins')
                   .when(hl.is_deletion(ht.alleles[0], ht.alleles[1]), 'del')
                   .default('complex')
                   )
    ht = ht.annotate(allele_data=ht.allele_data.annotate(allele_type=allele_type,
                                                         was_mixed=ht.allele_data.variant_type == 'mixed'))
    return ht


def generate_ac(mt: hl.MatrixTable, fam_file: str) -> hl.Table:
    """
    Creates Table with QC samples, QC samples removing children and release samples raw and adj ACs.
    """
    # mt = mt.filter_cols(mt.meta.high_quality)
    fam_ht = hl.import_fam(fam_file, delimiter="\t")
    mt = mt.annotate_cols(unrelated_sample=hl.is_missing(fam_ht[mt.s]))
    mt = mt.filter_rows(hl.len(mt.alleles) > 1)
    mt = mt.annotate_rows(
        ac_qc_samples_raw=hl.agg.sum(mt.GT.n_alt_alleles()),
        ac_qc_samples_adj=hl.agg.filter(
            mt.adj, hl.agg.sum(mt.GT.n_alt_alleles())),
    )
    return mt.rows()


def create_inbreeding_ht_with_ac_and_allele_data(varqc_mtfile: str, pedfile: str, inbreeding_htfile: str, qc_ac_htfile: str, allele_data_htfile: str):
    '''
    Inbreeding, allele data and qc_ac annotations
    :param str pedfile: Plink-format ped file
    :param str varqc_mtfile: Input matrixtable
    :param str inbreeding_htfile: Inbreeding hail table file
    :param str qc_ac_htfile: qc_ac hail table file
    :param str allele_data_htfile: Allele data htfile
    '''
    mt = hl.read_matrix_table(varqc_mtfile)
    # inbreeding ht
    mt_inbreeding = mt.annotate_rows(InbreedingCoeff=bi_allelic_site_inbreeding_expr(mt.GT))
    #mt = mt.key_rows_by('locus').distinct_by_row().key_rows_by('locus', 'alleles')
    ht_inbreeding = mt_inbreeding.rows()
    # allele data and qc_ac ht
    allele_data_ht = generate_allele_data(mt)
    qc_ac_ht = generate_ac(mt, pedfile)
    # write to file
    ht_inbreeding.write(inbreeding_htfile , overwrite=True)
    qc_ac_ht.write(qc_ac_htfile , overwrite=True)
    allele_data_ht.write(allele_data_htfile, overwrite=True)


def main():
    #set up
    args = get_options()
    inputs = parse_config()
    mtdir = inputs['matrixtables_lustre_dir']
    resourcedir = inputs['resource_dir']
    training_sets_dir = inputs['training_set_dir']

    # initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    # get truth set ht
    truth_ht_file = resourcedir + "truthset_table.ht"
    if args.truth or args.all:
        omni = training_sets_dir + "1000G_omni2.5.hg38.ht"
        mills = training_sets_dir + "Mills_and_1000G_gold_standard.indels.hg38.ht"
        thousand_genomes = training_sets_dir + "1000G_phase1.snps.high_confidence.hg38.ht"
        hapmap = training_sets_dir + "hapmap_3.3.hg38.ht"
        get_truth_ht(omni, mills, thousand_genomes, hapmap, truth_ht_file)

    #add hail variant QC
    if args.annotation or args.all:
        mtfile = mtdir + "mt_pops_QC_filters_sample_qc.shrank.mt"
        varqc_mtfile = mtdir + "mt_varqc.mt"
        varqc_mtfile_split = mtdir + "mt_varqc_splitmulti.mt"

        split_multi_and_var_qc(mtfile, varqc_mtfile, varqc_mtfile_split)
        pedfile = 'file:///lustre/scratch123/projects/gnh_industry/Genes_and_Health_2023_02_44k/GH_44k_668-trios_QCed.mercury.consistent.fam'

        #get complete trios, family annotation, dnm annotation
        trio_mtfile = mtdir + "trios.mt"
        trio_stats_htfile = mtdir + "trio_stats.ht"
        fam_stats_htfile = mtdir + "family_stats.ht"
        fam_stats_mtfile = mtdir + "family_stats.mt"
        fam_stats_gnomad_mtfile = mtdir + "family_stats_gnomad.mt"
        gnomad_htfile = resourcedir + "gnomad.exomes.r2.1.1.sites.liftover_grch38.ht"
        dnm_htfile = mtdir + "denovo_table.ht"
        trio_family_dnm_annotation(varqc_mtfile_split, pedfile, trio_mtfile, trio_stats_htfile, fam_stats_htfile, fam_stats_mtfile, fam_stats_gnomad_mtfile, gnomad_htfile, dnm_htfile)

        #create inbreeding ht
        inbreeding_htfile = mtdir + "inbreeding.ht"
        qc_ac_htfile = mtdir + "qc_ac.ht"
        allele_data_htfile = mtdir + "allele_data.ht"
        create_inbreeding_ht_with_ac_and_allele_data(varqc_mtfile, pedfile, inbreeding_htfile, qc_ac_htfile, allele_data_htfile)
        rm_mt(varqc_mtfile)


if __name__ == '__main__':
    main()
