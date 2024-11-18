"""
All fucntions used to filter samples/variants
"""

import hail as hl


def filter_matrix_for_ldprune(mt: hl.MatrixTable, long_range_ld_file: str) -> hl.MatrixTable:
    """
    This fucntion filters samples and variants to perform LD pruning,
    pruning of related samples

    The fucntion is taken from the step 2/3-annotate_and_filter
    In future step 2/3 should be refactored to use this function
    """
    # use only autosomes
    mt = mt.filter_rows(mt.locus.in_autosome())
    # split multiallelic variants and remove them
    mt = hl.split_multi_hts(mt)  # this shouldn't do anything as only biallelic sites are used
    mt = mt.filter_rows(mt.was_split == True, keep=False)
    # keep only SNPs
    mt = mt.filter_rows(hl.is_snp(mt.alleles[0], mt.alleles[1]))
    # keep good variants using hail variant_qc and thre filters
    mt_vqc = hl.variant_qc(mt, name="variant_QC_Hail")
    mt_vqc_filtered = mt_vqc.filter_rows(
        (mt_vqc.variant_QC_Hail.call_rate >= 0.99)
        & (mt_vqc.variant_QC_Hail.AF[1] >= 0.05)
        & (mt_vqc.variant_QC_Hail.p_value_hwe >= 1e-5)
    )
    # remove long ld regions
    long_range_ld_to_exclude = hl.import_bed(long_range_ld_file, reference_genome="GRCh38")
    mt_vqc_filtered = mt_vqc_filtered.filter_rows(
        hl.is_defined(long_range_ld_to_exclude[mt_vqc_filtered.locus]), keep=False
    )
    # remove palindromes
    mt_non_pal = mt_vqc_filtered.filter_rows(
        (mt_vqc_filtered.alleles[0] == "G") & (mt_vqc_filtered.alleles[1] == "C"), keep=False
    )
    mt_non_pal = mt_non_pal.filter_rows((mt_non_pal.alleles[0] == "C") & (mt_non_pal.alleles[1] == "G"), keep=False)
    mt_non_pal = mt_non_pal.filter_rows((mt_non_pal.alleles[0] == "A") & (mt_non_pal.alleles[1] == "T"), keep=False)
    mt_non_pal = mt_non_pal.filter_rows((mt_non_pal.alleles[0] == "T") & (mt_non_pal.alleles[1] == "A"), keep=False)

    return mt_non_pal
