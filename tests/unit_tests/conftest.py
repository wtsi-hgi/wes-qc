import pytest
import hail as hl


@pytest.fixture(scope="session")
def mt_sex_chr():
    mt = hl.balding_nichols_model(n_populations=1, n_samples=1, n_variants=50, reference_genome="GRCh38")
    x_mt = hl.balding_nichols_model(n_populations=1, n_samples=1, n_variants=25, reference_genome="GRCh38")
    return mt.union_rows(x_mt)


@pytest.fixture(scope="session")
def variation_mt(n_vars=32, include_indels=True, include_multiallelic=True):
    snvs = [["A", "T"], ["G", "C"]]
    indels = [["A", "ATA"], ["ATAT", "A"]]
    snvs_multallele = [["A", "T", "G"], ["G", "C", "T"]]
    indels_multallele = [["A", "ATAT", "AT"], ["ATAT", "A", "AT"]]

    variations = snvs
    if include_indels:
        variations.extend(indels)
    if include_multiallelic:
        variations.extend(snvs_multallele)
    if include_indels and include_multiallelic:
        variations.extend(indels_multallele)

    variations_all = variations * (n_vars // len(variations) + 1)

    alleles_list = [hl.literal(change) for change in variations_all[:n_vars]]
    loci = [
        hl.locus_from_global_position(hl.rand_int64(1, hl.get_reference("GRCh38").lengths["chr1"]), "GRCh38")
        for i in range(len(alleles_list))
    ]

    rows = hl.Table.parallelize([hl.struct(locus=locus, alleles=alleles) for locus, alleles in zip(loci, alleles_list)])

    mt = hl.MatrixTable.from_rows_table(
        rows,
    ).key_rows_by("locus", "alleles")

    mt = mt.annotate_cols(s="Mock_sample")
    mt = mt.annotate_entries(GT=hl.call(hl.rand_bool(0.5), hl.rand_bool(0.5)))

    return mt


"""
def generate_indel_mt():
    mt = hl.balding_nichols_model(n_populations=3, n_samples=1, n_variants=10, reference_genome="GRCh38")
    mt = mt.drop("ancestral_af", "af")
    loci = [
        hl.locus_from_global_position(hl.rand_int64(1, hl.get_reference("GRCh38").lengths["chr1"]), "GRCh38")
        for i in range(100)
    ]
    alleles = [hl.literal(["A", "ATAT", "AT"]), hl.literal(["ATAT", "A", "AT"])]
    indels = hl.Table.parallelize([hl.struct(locus=locus, alleles=random.choice(alleles)) for locus in loci])
    indels = indels.key_by("locus", "alleles")
    return mt.rows().union(indels)
"""
