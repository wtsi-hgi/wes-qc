"""
All fucntions used to filter samples/variants
"""

import hail as hl  # type: ignore


def filter_matrix_for_ldprune(
    mt: hl.MatrixTable,
    long_range_ld_file: str,
    call_rate_threshold: float = 0.99,
    af_threshold: float = 0.05,
    hwe_threshold: float = 1e-5,
) -> hl.MatrixTable:
    """
    This function filters samples and variants to perform LD pruning,
    pruning of related samples

    The function is taken from the step 2/3-annotate_and_filter
    In future step 2/3 should be refactored to use this function
    """
    mt = filter_mt_autosome_biallelic_snvs(mt)
    mt_vqc_filtered = filter_vars_for_quality(mt, af_threshold, call_rate_threshold, hwe_threshold)
    mt_vqc_filtered = remove_ld_regions(mt_vqc_filtered, long_range_ld_file)
    mt_non_pal = remove_palindromes(mt_vqc_filtered)
    return mt_non_pal


def filter_mt_autosome_biallelic_snvs(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Keeps only bi-allelic autosome SNVs
    """
    mt = mt.filter_rows(mt.locus.in_autosome())
    # split multiallelic variants and remove them
    mt = hl.split_multi_hts(mt)  # this shouldn't do anything as only biallelic sites are used
    mt = mt.filter_rows(mt.was_split is True, keep=False)
    # keep only SNVs
    mt = mt.filter_rows(hl.is_snp(mt.alleles[0], mt.alleles[1]))
    return mt


def filter_vars_for_quality(
    mt: hl.MatrixTable, af_threshold: float, call_rate_threshold: float, hwe_threshold: float
) -> hl.MatrixTable:
    """
    Keeps "good" variants passing the thresholds for QC metrics
    """
    print(f"\nDEBUG === Total variants: {mt.count_rows()}\n")
    mt_vqc = hl.variant_qc(mt, name="variant_QC_Hail")
    mt_vqc_filtered = mt_vqc.filter_rows(
        (mt_vqc.variant_QC_Hail.call_rate >= call_rate_threshold)
        & (mt_vqc.variant_QC_Hail.AF[1] >= af_threshold)
        & (mt_vqc.variant_QC_Hail.p_value_hwe >= hwe_threshold)
    )
    print(f"\nDEBUG === Total variants survived after quality filtering: {mt_vqc_filtered.count_rows()}\n")
    return mt_vqc_filtered


def remove_ld_regions(mt: hl.MatrixTable, long_range_ld_file: str) -> hl.MatrixTable:
    # remove long ld regions
    long_range_ld_to_exclude = hl.import_bed(long_range_ld_file, reference_genome="GRCh38")
    mt = mt.filter_rows(hl.is_defined(long_range_ld_to_exclude[mt.locus]), keep=False)
    return mt


def remove_palindromes(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Removes palindromic SNVs from the matrix table
    A palindromic SNP (also known as an "ambiguous" SNP) is a SNP
    in which the possible alleles for that SNP are the same alleles
    that would pair with each other in the double helix structure.
    e.g., C/G on forward is G/C on the reverse
    https://mr-dictionary.mrcieu.ac.uk/term/palindrome/

    To get rid of these SNVs we remove all paired combinations: A<->T, G<->C
    """
    mt_non_pal = mt.filter_rows((mt.alleles[0] == "G") & (mt.alleles[1] == "C"), keep=False)
    mt_non_pal = mt_non_pal.filter_rows((mt_non_pal.alleles[0] == "C") & (mt_non_pal.alleles[1] == "G"), keep=False)
    mt_non_pal = mt_non_pal.filter_rows((mt_non_pal.alleles[0] == "A") & (mt_non_pal.alleles[1] == "T"), keep=False)
    mt_non_pal = mt_non_pal.filter_rows((mt_non_pal.alleles[0] == "T") & (mt_non_pal.alleles[1] == "A"), keep=False)
    return mt_non_pal
