#apply a range of different hard filters (RF bin and genotype) to SNPs and indels, also add gene and consequence
#removes samples which fail identity checks
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
    args = parser.parse_args()
    if not args.runhash:
        print("--runhash must be specified")
        exit(1)

    return args


def remove_samples(mt: hl.MatrixTable, exclude_file: str):
    '''
    Remove samples in file of samples which fail identity checks
    :param hl.MatrixTable: Input MatrixTable
    :param str exclude_file: path of file with samples to exclude
    :return: hl.MatrixTable
    '''
    excl_ht = hl.import_table(exclude_file, impute=True, key = 'exomeID')
    mt = mt.filter_cols(hl.is_defined(excl_ht[mt.s]), keep=False)

    return mt


def annotate_cq_rf(mt: hl.MatrixTable, rf_htfile: str, cqfile: str) -> hl.MatrixTable:
    '''
    Annotate with RF bin, consequence, gene name, hgnc id, and pass/fail for 3 combinations of filters for SNPs,
    3 combinations of filters for indels and missingness for each set of filters
    :param hl.MatrixTable mtfile: Input matrixtable
    :param str rf_htfile: random forest hail table file
    :param str cqfile: consequence file
    :return hl.matrixTable:
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

    cq_ht = hl.import_table(cqfile, types={'f0': 'str', 'f1': 'int32', 'f2': 'str',
                            'f3': 'str', 'f4': 'str', 'f5': 'str', 'f6': 'str', 'f7': 'str', 'f8': 'str'}, no_header=True)
    cq_ht = cq_ht.annotate(chr=cq_ht.f0)
    cq_ht = cq_ht.annotate(pos=cq_ht.f1)
    cq_ht = cq_ht.annotate(rs=cq_ht.f2)
    cq_ht = cq_ht.annotate(ref=cq_ht.f3)
    cq_ht = cq_ht.annotate(alt=cq_ht.f4)
    cq_ht = cq_ht.annotate(CSQ=cq_ht.f5)
    cq_ht = cq_ht.annotate(consequence=cq_ht.f6)
    cq_ht = cq_ht.annotate(gene=cq_ht.f7)
    cq_ht = cq_ht.annotate(hgnc_id=cq_ht.f8)
    cq_ht = cq_ht.key_by(
        locus=hl.locus(cq_ht.chr, cq_ht.pos), alleles=[cq_ht.ref, cq_ht.alt])
    cq_ht = cq_ht.drop(cq_ht.f0, cq_ht.f1, cq_ht.f2, cq_ht.f3, cq_ht.f4, cq_ht.chr, cq_ht.pos, cq_ht.ref, cq_ht.alt)
    cq_ht = cq_ht.key_by(cq_ht.locus, cq_ht.alleles)
    mt = mt.annotate_rows(
        info=mt.info.annotate(
            CSQ=cq_ht[mt.row_key].CSQ,
            consequence=cq_ht[mt.row_key].consequence,
            gene=cq_ht[mt.row_key].gene,
            hgnc_id=cq_ht[mt.row_key].hgnc_id)
    )

    return mt


def apply_hard_filters(mt: hl.MatrixTable, hard_filters: dict) -> hl.MatrixTable:
    '''
    Apply hard filters and annotate missingness
    :param hl.MatrixTable mt: Input MatrixTable
    :param dict hard_filters: Filters to apply
    :return hl.MatrixTable:
    '''

    # stringent_condition = ( (
    #                             (hl.is_snp(mt.alleles[0], mt.alleles[1])) 
    #                             & (mt.info.rf_bin <= hard_filters['snp']['stringent']['bin'])
    #                             & (mt.DP < hard_filters['snp']['stringent']['dp']) 
    #                             & (mt.GQ < hard_filters['snp']['stringent']['gq']) 
    #                             & ( (mt.GT.is_het() & (mt.HetAB < hard_filters['snp']['stringent']['ab'])) )
    #                             ) |
    #                             (
    #                             (hl.is_indel(mt.alleles[0], mt.alleles[1])) 
    #                             & (mt.info.rf_bin <= hard_filters['indel']['stringent']['bin'])
    #                             & (mt.DP < hard_filters['indel']['stringent']['dp']) 
    #                             & (mt.GQ < hard_filters['indel']['stringent']['gq']) 
    #                             & ( (mt.GT.is_het() & (mt.HetAB < hard_filters['indel']['stringent']['ab'])) )
    #                             )   
    #                         )


    # medium_condition = ( ((hl.is_snp(mt.alleles[0], mt.alleles[1])) 
    #                             & (mt.info.rf_bin <= hard_filters['snp']['medium']['bin'])
    #                             & (mt.DP < hard_filters['snp']['medium']['dp']) 
    #                             & (mt.GQ < hard_filters['snp']['medium']['gq']) 
    #                             & ( (mt.GT.is_het() & (mt.HetAB < hard_filters['snp']['medium']['ab'])) ) 
    #                             ) | 
    #                             (
    #                             (hl.is_indel(mt.alleles[0], mt.alleles[1])) 
    #                             & (mt.info.rf_bin <= hard_filters['indel']['medium']['bin'])
    #                             & (mt.DP < hard_filters['indel']['medium']['dp']) 
    #                             & (mt.GQ < hard_filters['indel']['medium']['gq']) 
    #                             & ( (mt.GT.is_het() & (mt.HetAB < hard_filters['indel']['medium']['ab'])) )
    #                             )
    #                         )
    
    # relaxed_condition = ( ((hl.is_snp(mt.alleles[0], mt.alleles[1])) 
    #                             & (mt.info.rf_bin <= hard_filters['snp']['relaxed']['bin'])
    #                             & (mt.DP < hard_filters['snp']['relaxed']['dp']) 
    #                             & (mt.GQ < hard_filters['snp']['relaxed']['gq']) 
    #                             & ( (mt.GT.is_het() & (mt.HetAB < hard_filters['snp']['relaxed']['ab'])) )
    #                             ) |
    #                             (
    #                             (hl.is_indel(mt.alleles[0], mt.alleles[1])) 
    #                             & (mt.info.rf_bin <= hard_filters['indel']['relaxed']['bin'])
    #                             & (mt.DP < hard_filters['indel']['relaxed']['dp']) 
    #                             & (mt.GQ < hard_filters['indel']['relaxed']['gq']) 
    #                             & ( (mt.GT.is_het() & (mt.HetAB < hard_filters['indel']['relaxed']['ab'])) ) 
    #                             )                               
    #                         )

    # mt = mt.annotate_entries(
    #     stringent_filters = hl.if_else(stringent_condition, 'Fail', 'Pass'),
    #     medium_filters = hl.if_else(medium_condition, 'Fail', 'Pass'),
    #     relaxed_filters = hl.if_else(relaxed_condition, 'Fail', 'Pass')
    # )

    stringent_condition = (
        (
        (hl.is_snp(mt.alleles[0], mt.alleles[1]))
        & (mt.info.rf_bin <= hard_filters['snp']['stringent']['bin'])
        & (mt.DP >= hard_filters['snp']['stringent']['dp'])
        & (mt.GQ >= hard_filters['snp']['stringent']['gq'])
        & (
            (mt.GT.is_het() & (mt.HetAB >= hard_filters['snp']['stringent']['ab'])) |
            (mt.GT.is_hom_ref()) |
            (mt.GT.is_hom_var())
        )

    ) |
        (
        (hl.is_indel(mt.alleles[0], mt.alleles[1]))
        & (mt.info.rf_bin <= hard_filters['indel']['stringent']['bin'])
        & (mt.DP >= hard_filters['indel']['stringent']['dp'])
        & (mt.GQ >= hard_filters['indel']['stringent']['gq'])
        & (
            (mt.GT.is_het() & (mt.HetAB >= hard_filters['indel']['stringent']['ab'])) |
            (mt.GT.is_hom_ref()) |
            (mt.GT.is_hom_var())
        )

    )
    )


    medium_condition = (
        (
        (hl.is_snp(mt.alleles[0], mt.alleles[1]))
        & (mt.info.rf_bin <= hard_filters['snp']['medium']['bin'])
        & (mt.DP >= hard_filters['snp']['medium']['dp'])
        & (mt.GQ >= hard_filters['snp']['medium']['gq'])
        & (
            (mt.GT.is_het() & (mt.HetAB >= hard_filters['snp']['medium']['ab'])) |
            (mt.GT.is_hom_ref()) |
            (mt.GT.is_hom_var())
        )

    ) |
        (
        (hl.is_indel(mt.alleles[0], mt.alleles[1]))
        & (mt.info.rf_bin <= hard_filters['indel']['medium']['bin'])
        & (mt.DP >= hard_filters['indel']['medium']['dp'])
        & (mt.GQ >= hard_filters['indel']['medium']['gq'])
        & (
            (mt.GT.is_het() & (mt.HetAB >= hard_filters['indel']['medium']['ab'])) |
            (mt.GT.is_hom_ref()) |
            (mt.GT.is_hom_var())
        )

    )
    )
    
    relaxed_condition = (
        (
        (hl.is_snp(mt.alleles[0], mt.alleles[1]))
        & (mt.info.rf_bin <= hard_filters['snp']['relaxed']['bin'])
        & (mt.DP >= hard_filters['snp']['relaxed']['dp'])
        & (mt.GQ >= hard_filters['snp']['relaxed']['gq'])
        & (
            (mt.GT.is_het() & (mt.HetAB >= hard_filters['snp']['relaxed']['ab'])) |
            (mt.GT.is_hom_ref()) |
            (mt.GT.is_hom_var())
        )

    ) |
        (
        (hl.is_indel(mt.alleles[0], mt.alleles[1]))
        & (mt.info.rf_bin <= hard_filters['indel']['relaxed']['bin'])
        & (mt.DP >= hard_filters['indel']['relaxed']['dp'])
        & (mt.GQ >= hard_filters['indel']['relaxed']['gq'])
        & (
            (mt.GT.is_het() & (mt.HetAB >= hard_filters['indel']['relaxed']['ab'])) |
            (mt.GT.is_hom_ref()) |
            (mt.GT.is_hom_var())
        )

    )
    )

    mt = mt.annotate_entries(
        stringent_filters = hl.if_else(stringent_condition, 'Pass', 'Fail'),
        medium_filters = hl.if_else(medium_condition, 'Pass', 'Fail'),
        relaxed_filters = hl.if_else(relaxed_condition, 'Pass', 'Fail')
    )
    #annotate variants with fraction passing/failing each set of filters
    mt = mt.annotate_rows(stringent_pass_count = hl.agg.count_where(mt.stringent_filters == 'Pass'))
    mt = mt.annotate_rows(stringent_fail_count = hl.agg.count_where(mt.stringent_filters == 'Fail'))
    mt = mt.annotate_rows(info=mt.info.annotate(fraction_pass_stringent_filters = mt.stringent_pass_count/(mt.stringent_pass_count + mt.stringent_fail_count)))
    
    mt = mt.annotate_rows(medium_pass_count = hl.agg.count_where(mt.medium_filters == 'Pass'))
    mt = mt.annotate_rows(medium_fail_count = hl.agg.count_where(mt.medium_filters == 'Fail'))
    mt = mt.annotate_rows(info=mt.info.annotate(fraction_pass_medium_filters = mt.medium_pass_count/(mt.medium_pass_count + mt.medium_fail_count)))
    
    mt = mt.annotate_rows(relaxed_pass_count = hl.agg.count_where(mt.relaxed_filters == 'Pass'))
    mt = mt.annotate_rows(relaxed_fail_count = hl.agg.count_where(mt.relaxed_filters == 'Fail'))
    mt = mt.annotate_rows(info=mt.info.annotate(fraction_pass_relaxed_filters = mt.relaxed_pass_count/(mt.relaxed_pass_count + mt.relaxed_fail_count)))

    return mt


def main():
    #set up
    args = get_options()
    inputs = parse_config()
    mtdir = inputs['matrixtables_lustre_dir']
    rf_dir = inputs['var_qc_rf_dir']
    resourcedir = inputs['resource_dir']
    hard_filters = inputs['hard_filters']
    annotdir = inputs['annotation_lustre_dir']

    # initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    exclude_file = annotdir + "to_be_excluded_exome.txt"
    rf_htfile = rf_dir + args.runhash + "/_gnomad_score_binning_tmp.ht"
    mtfile = mtdir + "mt_varqc_splitmulti.mt"
    cqfile = resourcedir + "all_consequences_with_gene_and_csq.txt"
    mtfile_annot = mtdir + "mt_hard_filter_combinations.mt"

    mt = hl.read_matrix_table(mtfile)

    #remove unwanted samples
    mt = remove_samples(mt, exclude_file)

    #annotate mt with consequence, gene, rf bin
    mt_annot = annotate_cq_rf(mt, rf_htfile, cqfile)

    #annotate with all combinations of filters (pass/fail) and add missingness
    mt_annot = apply_hard_filters(mt_annot, hard_filters)
    mt_annot.write(mtfile_annot, overwrite = True)


if __name__ == '__main__':
    main()