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


@pytest.fixture(scope="session")
def pca_scores_table() -> hl.Table:
    samples = [
        hl.Struct(
            s="NA20894",
            scores=[-0.06455414043179736, -0.14946688518899134, -0.2456780547093593, -0.23921033922105084],
            known_pop="SAS",
        ),
        hl.Struct(
            s="NA20903",
            scores=[-0.07126483958325104, -0.24364068671180367, -0.16253488146259326, -0.27010859769162615],
            known_pop=None,
        ),
        hl.Struct(
            s="NA19751",
            scores=[-0.15598650919318002, -0.0442263377892484, 0.13368269478712813, -0.04503040877066483],
            known_pop="AMR",
        ),
        hl.Struct(
            s="NA20764",
            scores=[-0.011736033892254222, -0.34824623124705745, -0.04832096176489377, -0.08749133180605716],
            known_pop="EUR",
        ),
        hl.Struct(
            s="NA20766",
            scores=[-0.0944388589031874, -0.2037407129001564, -0.0036595171865637217, 0.07110249529033492],
            known_pop="EUR",
        ),
        hl.Struct(
            s="NA20772",
            scores=[-0.1031612294892679, -0.29375836490379226, -0.00945596551086085, -0.025279767480578638],
            known_pop=None,
        ),
        hl.Struct(
            s="NA19625",
            scores=[0.18882141015602674, 0.04675680416694248, 0.12086001948711221, 0.05616217397859573],
            known_pop="AFR",
        ),
        hl.Struct(
            s="NA19082",
            scores=[-0.10142021239596818, 0.2640727135985397, -0.156534646894955, 0.21255105227110127],
            known_pop="EAS",
        ),
    ]

    # Convert the list to a Hail table
    return hl.Table.parallelize(samples)


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
