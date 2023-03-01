#
import hail as hl
import pyspark
import argparse
from utils.utils import parse_config


def get_options():
    '''
    Get options from the command line
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("--runhash", help="RF training run hash")
    args = parser.parse_args()
    if not args.runhash:
        print("--runhash must be specified")
        exit(1)

    return args


def add_cq_annotation(htfile: str, synonymous_file: str, ht_cq_file: str):
    '''
    Add annotation for synonymous variants
    ;param str htfile: Random forest hail table file
    ;param str synonymous_file: text file containing synonymous variants
    ;param str ht_cq_fle: Output hail table file annotated with synonymous variants
    '''
    synonymous_ht=hl.import_table(synonymous_file,types={'f0':'str','f1':'int32', 'f2':'str','f3':'str','f4':'str', 'f5':'str'}, no_header=True)
    synonymous_ht=synonymous_ht.annotate(chr=synonymous_ht.f0)
    synonymous_ht=synonymous_ht.annotate(pos=synonymous_ht.f1)
    synonymous_ht=synonymous_ht.annotate(rs=synonymous_ht.f2)
    synonymous_ht=synonymous_ht.annotate(ref=synonymous_ht.f3)
    synonymous_ht=synonymous_ht.annotate(alt=synonymous_ht.f4)
    synonymous_ht=synonymous_ht.annotate(consequence=synonymous_ht.f5)
    synonymous_ht = synonymous_ht.key_by(
        locus=hl.locus(synonymous_ht.chr, synonymous_ht.pos), alleles=[synonymous_ht.ref,synonymous_ht.alt])
    synonymous_ht=synonymous_ht.drop(synonymous_ht.f0,synonymous_ht.f1,synonymous_ht.f2,synonymous_ht.f3,synonymous_ht.f4,synonymous_ht.chr,synonymous_ht.pos,synonymous_ht.ref,synonymous_ht.alt)
    synonymous_ht = synonymous_ht.key_by(synonymous_ht.locus, synonymous_ht.alleles)

    ht = hl.read_table(htfile)
    ht=ht.annotate(consequence=synonymous_ht[ht.key].consequence)
    ht.write(ht_cq_file, overwrite=True)


def dnm_and_family_annotation(ht_cq_fle: str, dnm_htfile: str, fam_stats_htfile: str, trio_stats_htfile: str, family_annot_htfile: str):
    '''
    Annotate RF result hail table with DNMs and family stats
    ;param str ht_cq_fle: Input hail table file annotated with synonymous variants
    :param str dnm_htfile: De novo hail table file
    :param str fam_stats_htfile: Family stats hail table file
    :param str trio_stats_htfile: Trio stats hail table file
    :param str family_annot_htfile: Output hail table annotated with DNMs and family stats file
    '''
    ht = hl.read_table(ht_cq_fle)
    dnm_ht = hl.read_table(dnm_htfile)
    fam_stats_ht = hl.read_table(fam_stats_htfile)
    trio_stats_ht = hl.read_table(trio_stats_htfile)
    ht=ht.annotate(de_novo_data=dnm_ht[ht.key].de_novo_data)
    ht=ht.annotate(family_stats=fam_stats_ht[ht.key].family_stats)
    ht=ht.annotate(fam=trio_stats_ht[ht.key])
    ht.write(family_annot_htfile, overwrite=True)


def count_trans_untransmitted_singletons(mt_filtered: hl.MatrixTable, ht: hl.Table) -> hl.Table:
    '''
    Count untransmitted and transmitted singletons
    :param hl.MatrixTable mt_filtered: Filtered trio matrixtable
    :param hl.Table ht: Output hail table
    '''
    

    # mt_trans = mt_filtered.filter_entries(mt_filtered.variant_qc.AC[1] == 2)
    # mt_untrans = mt_filtered.filter_entries(mt_filtered.variant_qc.AC[1] == 1)
    mt_trans = mt_filtered.filter_entries(mt_filtered.varqc_trios.AC[1] == 2)
    mt_untrans = mt_filtered.filter_entries(mt_filtered.varqc_trios.AC[1] == 1)
    
    mt_trans_count=mt_trans.group_cols_by(mt_trans.id).aggregate(transmitted_singletons_count=hl.agg.count_where(
                               # (mt_trans.info.AC[0] == 2) &
                                (mt_trans.proband_entry.GT.is_non_ref()) &
                                (
                                (mt_trans.father_entry.GT.is_non_ref()  )
                                 |
                                (mt_trans.mother_entry.GT.is_non_ref())
                                )
                                ))
    

    
    Total_transmitted_singletons=mt_trans_count.aggregate_entries(hl.agg.count_where(mt_trans_count.transmitted_singletons_count >0))
    print(Total_transmitted_singletons)
    
    mt_untrans_count = (mt_untrans.group_cols_by(mt_untrans.id).aggregate(
    untransmitted_singletons_count=hl.agg.count_where(
                   # (mt_untrans.info.AC[0] == 1) &
                    (mt_untrans.proband_entry.GT.is_hom_ref()) &
                    (
                    (mt_untrans.father_entry.GT.is_non_ref()) 
                    |
                    (mt_untrans.mother_entry.GT.is_non_ref())
                    )
                     )))
    Total_untransmitted_singletons=mt_untrans_count.aggregate_entries(hl.agg.count_where(mt_untrans_count.untransmitted_singletons_count >0))
    print(Total_untransmitted_singletons)
    
    print(f"\nTransmitted singletons:{Total_transmitted_singletons}\n")
    print(f"\nUntransmitted singletons:{Total_untransmitted_singletons}")

    Ratio_transmitted_untransmitted=Total_transmitted_singletons/Total_untransmitted_singletons
    print(Ratio_transmitted_untransmitted)
    print(f"\nRatio:{Ratio_transmitted_untransmitted}\n")
    mt2=mt_trans_count.annotate_rows(variant_transmitted_singletons=hl.agg.count_where(mt_trans_count.transmitted_singletons_count==1))
    mt2.variant_transmitted_singletons.summarize()

    mt3=mt_untrans_count.annotate_rows(variant_untransmitted_singletons=hl.agg.count_where(mt_untrans_count.untransmitted_singletons_count==1))
    mt3.variant_untransmitted_singletons.summarize()

    ht=ht.annotate(variant_transmitted_singletons=mt2.rows()[ht.key].variant_transmitted_singletons)
    ht=ht.annotate(variant_untransmitted_singletons=mt3.rows()[ht.key].variant_untransmitted_singletons)
    return(ht)    


def transmitted_singleton_annotation(family_annot_htfile: str, trio_mtfile: str, trio_filtered_mtfile: str, trans_sing_htfile: str):
    '''
    Annotate MT with transmited singletons and create final variant QC HT for ranking
    :param str family_annot_htfile: Family annotation hail table file
    :param str trio_mtfile: Trio annotation hail matrixtable file
    :param str trio_filtered_mtfile: Trio filtered hail matrixtable file
    :param str trans_sing_htfile: Variant QC hail table file with transmitted singleton annotation
    '''
    print("Annotating with transmitted singletons")
    ht = hl.read_table(family_annot_htfile)
    mt_trios = hl.read_matrix_table(trio_mtfile)
    mt_trios = mt_trios.annotate_rows(consequence=ht[mt_trios.row_key].consequence)
    #there is a save step here in Pavlos file, is it needed? 
    # mt_filtered = mt_trios.filter_rows((mt_trios.variant_qc.AC[1] <= 2) & (mt_trios.consequence == "synonymous_variant"))
    # mt_filtered = mt_trios.filter_entries((mt_trios.variant_qc.AC[1] <= 2) & (mt_trios.consequence == "synonymous_variant"))
    mt_filtered = mt_trios.filter_rows((mt_trios.varqc_trios.AC[1] <= 2) & (mt_trios.consequence == "synonymous_variant"))
    mt_filtered = mt_trios.filter_entries((mt_trios.varqc_trios.AC[1] <= 2) & (mt_trios.consequence == "synonymous_variant"))
    mt_filtered.write(trio_filtered_mtfile, overwrite=True)

    ht = count_trans_untransmitted_singletons(mt_filtered, ht)
    ht.write(trans_sing_htfile, overwrite=True)    


# def run_tdt(mtfile: str, trans_sing_htfile: str, pedfile: str, tdt_htfile: str):
#     '''
#     Run transmission disequilibrium test to get counts of number transmitted and number unstransmitted for each 
#     variant and annotate the RF output file with this
#     :param str mtfile: Mtfile used to create random forest input table
#     :param str trans_sing_htfile: Variant QC hail table file with transmitted singleton annotation
#     :param str pedfile: Pedfile
#     :param str tdt_htfile: Htfile with transmitted/untransmitted counts
#     '''
#     print("Annotating with transmitted/unstransmitted counts")
#     mt = hl.read_matrix_table(mtfile)
#     pedigree = hl.Pedigree.read(pedfile)
#     tdt_ht = hl.transmission_disequilibrium_test(mt, pedigree)
#     ht = hl.read_table(trans_sing_htfile)
#     ht=ht.annotate(n_transmitted=tdt_ht[ht.key].t)
#     ht=ht.annotate(n_untransmitted=tdt_ht[ht.key].u)
#     ht.write(tdt_htfile, overwrite = True)


def annotate_gnomad(tdt_htfile: str, gnomad_htfile: str, final_htfile: str):
    '''
    Annotate with gnomad allele frequencies
    :param str tdt_htfile: Htfile with transmitted/untransmitted counts
    :param str gnomad_htfile: Gnomad annotation hail table file
    :param str final_htfile: Final RF htfile for ranking and binning
    '''
    print("Annotating with gnomad AF")
    ht = hl.read_table(tdt_htfile)
    gnomad_ht = hl.read_table(gnomad_htfile)
    ht = ht.annotate(gnomad_af=gnomad_ht[ht.key].freq[0].AF)
    ht.write(final_htfile, overwrite = True)


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

    htfile = rf_dir + args.runhash + "/rf_result.ht"
    #annotate with synonymous CQs
    synonymous_file = resourcedir + "synonymous_variants.txt"
    ht_cq_file = rf_dir + args.runhash + "/rf_result_with_synonymous.ht"
    add_cq_annotation(htfile, synonymous_file, ht_cq_file)
    #annotate with family stats and DNMs
    dnm_htfile = mtdir + "denovo_table.ht"
    fam_stats_htfile = mtdir + "family_stats.ht"
    trio_stats_htfile = mtdir + "trio_stats.not-bi-only.ht"
    family_annot_htfile = rf_dir + args.runhash + "/rf_result_denovo_family_stats.ht"
    dnm_and_family_annotation(ht_cq_file, dnm_htfile, fam_stats_htfile, trio_stats_htfile, family_annot_htfile)

    #annotate with transmitted singletons
    trio_mtfile = mtdir + "trios.mt"
    trio_filtered_mtfile = mtdir + "trios_filtered_.mt"
    trans_sing_htfile = rf_dir + args.runhash + "/rf_result_trans_sing.ht"
    transmitted_singleton_annotation(family_annot_htfile, trio_mtfile, trio_filtered_mtfile, trans_sing_htfile)

    #annotate with number transmitted, number untransmitted (from transmission disequilibrium test)
    # mtfile = mtdir + "mt_varqc_splitmulti.mt"
    # pedfile = resourcedir + "trios.ped"
    # tdt_htfile = rf_dir + args.runhash + "/rf_result_tdt.ht"
    # run_tdt(mtfile, trans_sing_htfile, pedfile, tdt_htfile)

    #annotate with gnomad AF
    final_htfile = rf_dir + args.runhash + "/rf_result_final_for_ranking.ht"
    gnomad_htfile = resourcedir + "gnomad.exomes.r2.1.1.sites.liftover_grch38.ht"
    annotate_gnomad(trans_sing_htfile, gnomad_htfile, final_htfile)


if __name__ == '__main__':
    main()
