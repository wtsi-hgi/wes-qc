"""
All fucntions used to filter samples/variants
"""

import hail as hl


def filter_matrix_for_ldprune(
    mt: hl.MatrixTable,
    long_range_ld_file: str,
    call_rate_threshold: float = 0.99,
    af_threshold: float = 0.05,
    hwe_threshold: float = 1e-5,
) -> hl.MatrixTable:
    """
    This fucntion filters samples and variants to perform LD pruning,
    pruning of related samples

    The fucntion is taken from the step 2/3-annotate_and_filter
    In future step 2/3 should be refactored to use this function
    """
    mt = filter_mt_autosome_biallelic_snvs(mt)

    # keep good variants using hail variant_qc and thre filters
    mt_vqc = hl.variant_qc(mt, name="variant_QC_Hail")
    mt_vqc_filtered = mt_vqc.filter_rows(
        (mt_vqc.variant_QC_Hail.call_rate >= call_rate_threshold)
        & (mt_vqc.variant_QC_Hail.AF[1] >= af_threshold)
        & (mt_vqc.variant_QC_Hail.p_value_hwe >= hwe_threshold)
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


def filter_mt_autosome_biallelic_snvs(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    This function filters the matrix table for autosomes and SNVs
    """
    mt = mt.filter_rows(mt.locus.in_autosome())
    # split multiallelic variants and remove them
    mt = hl.split_multi_hts(mt)  # this shouldn't do anything as only biallelic sites are used
    mt = mt.filter_rows(mt.was_split is True, keep=False)
    # keep only SNVs
    mt = mt.filter_rows(hl.is_snp(mt.alleles[0], mt.alleles[1]))

    return mt
