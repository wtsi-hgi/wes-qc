import pytest
import hail as hl
import random


@pytest.fixture(scope="session")
def hail_data():
    mt = hl.balding_nichols_model(n_populations=3, n_samples=10, n_variants=50, reference_genome="GRCh38")
    x_mt = hl.balding_nichols_model(n_populations=3, n_samples=10, n_variants=50, reference_genome="GRCh38")
    indel_mt = generate_indel_mt()
    return mt.rows().union(x_mt.rows()).drop("ancestral_af", "af").union(indel_mt)


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
