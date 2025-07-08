import hail as hl  # type: ignore
from pytest import mark as m
import pytest
from wes_qc import filtering


def test_find_duplicated_variants_no_duplicates():
    # Create a sample MatrixTable with no duplicated variants
    rows_data = [
        {"locus": hl.Locus("1", 1000), "alleles": ["A", "T"]},
        {"locus": hl.Locus("2", 2000), "alleles": ["G", "C"]},
        {"locus": hl.Locus("3", 3000), "alleles": ["G", "A"]},
    ]

    rows = hl.Table.parallelize(rows_data, key=["locus", "alleles"])
    mt = hl.MatrixTable.from_rows_table(rows)

    # Call the function
    result = filtering.find_duplicated_variants(mt)

    # Verify the result is empty
    assert result.count() == 0


def test_find_duplicated_variants_multi_allelic_variants():
    # Create a sample MatrixTable with multi-allelic variants
    rows_data = [
        {"locus": hl.Locus("1", 1000), "alleles": ["A", "T", "G"]},
        {"locus": hl.Locus("1", 1000), "alleles": ["A", "T", "G"]},  # Duplicate
        {"locus": hl.Locus("1", 1000), "alleles": ["A", "T"]},  # Different alleles, not a duplicate
    ]

    rows = hl.Table.parallelize(rows_data, key=["locus", "alleles"])
    mt = hl.MatrixTable.from_rows_table(rows)
    result = filtering.find_duplicated_variants(mt)
    assert result.count() == 1
    row = result.collect()[0]
    assert row.variation[0] == hl.Locus("1", 1000)
    assert row.variation[1] == ["A", "T", "G"]
    assert row.var_count == 2


def test_find_duplicated_variants_multiple_duplicates():
    # Create a sample MatrixTable with multiple sets of duplicates
    rows_data = [
        {"locus": hl.Locus("1", 1000), "alleles": ["A", "T"]},
        {"locus": hl.Locus("1", 1000), "alleles": ["A", "T"]},  # Duplicate
        {"locus": hl.Locus("1", 1000), "alleles": ["A", "T"]},  # Duplicate
        {"locus": hl.Locus("2", 2000), "alleles": ["G", "C"]},
        {"locus": hl.Locus("2", 2000), "alleles": ["G", "C"]},  # Duplicate
        {"locus": hl.Locus("3", 3000), "alleles": ["G", "A"]},
    ]

    rows = hl.Table.parallelize(rows_data, key=["locus", "alleles"])
    mt = hl.MatrixTable.from_rows_table(rows)
    result = filtering.find_duplicated_variants(mt)
    assert result.count() == 2


@pytest.fixture
def sample_matrix_table():
    """Fixture to provide a Hail MatrixTable witj repeated samples."""
    n_samples = 10
    samples = hl.array([f"sample{i}" if i % 2 == 0 else "sample_dup" for i in range(n_samples)])
    mt = hl.utils.range_matrix_table(n_rows=5, n_cols=n_samples)
    mt = mt.annotate_cols(s=samples[mt.col_idx])
    return mt


def test_find_duplicated_samples_finds_duplicates(sample_matrix_table):
    """Test that the function identifies duplicated sample IDs."""
    duplicates = filtering.find_duplicated_samples(sample_matrix_table)
    assert isinstance(duplicates, hl.Table)
    assert duplicates.count() == 1  # we have only one sample with duplicated name
    col = duplicates.collect()[0]
    assert col.s == "sample_dup"
    assert col.s_count == 5


def test_find_duplicated_samples_no_duplicates():
    """Test when there are no duplicate sample IDs."""
    mt = hl.utils.range_matrix_table(n_rows=5, n_cols=5)
    mt = mt.annotate_cols(s=hl.str(mt.col_key))
    duplicates = filtering.find_duplicated_samples(mt)
    assert duplicates.count() == 0  # No duplicates expected


@m.context("When the input matrix table has variants on autosomes, SNVs, and is biallelic")
@m.it("should return a matrix table with only autosomes, SNVs, and non-split variants")
def test_filter_mt_autosome_biallelic_snvs(variation_mt, mt_sex_chr):
    # Apply the filter function
    filtered_mt = filtering.filter_mt_autosome_biallelic_snvs(mt_sex_chr)

    # Check that all variants are on autosomes
    assert filtered_mt.aggregate_rows(hl.agg.all(filtered_mt.locus.in_autosome()))

    filtered_mt = filtering.filter_mt_autosome_biallelic_snvs(variation_mt)

    # Check that all variants are SNVs
    assert filtered_mt.aggregate_rows(hl.agg.all(hl.is_snp(filtered_mt.alleles[0], filtered_mt.alleles[1])))

    # Check that we have no variants were split
    multiallelic_mt = filtered_mt.filter_rows(hl.len(filtered_mt.alleles) > 2)
    assert multiallelic_mt.count_rows() == 0


@m.context("Filters only variations with good quality metrics")
@m.it("Should run without exceptions")
def test_filter_vars_for_quality(mt_sex_chr):
    mt_vqc_filtered = filtering.filter_vars_for_quality(
        mt_sex_chr, call_rate_threshold=0.99, af_threshold=0.05, hwe_threshold=1e-5
    )
    mt_vqc_filtered.count_rows()


@pytest.fixture
def mock_hail_table():
    """Fixture to create a mock Hail Table"""
    data = [
        {"locus": hl.Locus("chr1", 1000, reference_genome="GRCh38")},
        {"locus": hl.Locus("chr1", 2000, reference_genome="GRCh38")},
        {"locus": hl.Locus("chr2", 3000, reference_genome="GRCh38")},
    ]
    ht = hl.Table.parallelize(data, hl.tstruct(locus=hl.tlocus(reference_genome="GRCh38")))
    ht = ht.key_by("locus")
    return ht


@pytest.fixture
def temp_bed(tmp_path):
    """Fixture to create a temporary BED file"""
    bed_content = "chr1\t500\t1500\nchr2\t2500\t3500\n"
    bed_file_name = tmp_path / "test.bed"
    bed_file_name.write_text(bed_content)
    bed = hl.import_bed(str(bed_file_name), reference_genome="GRCh38")
    bed = bed.key_by("interval")
    return bed


def test_filter_by_bed_filters_correctly(mock_hail_table, temp_bed):
    """Test that filter_by_bed correctly filters the table using the BED file"""
    result = filtering.filter_by_bed(mock_hail_table, temp_bed)
    expected_loci = [
        hl.Locus("chr1", 1000, reference_genome="GRCh38"),
        hl.Locus("chr2", 3000, reference_genome="GRCh38"),
    ]
    result_loci = result.locus.collect()
    assert result_loci == expected_loci


def test_filter_by_bed_no_matches(mock_hail_table, tmp_path):
    """Test that filter_by_bed handles cases where no loci match"""
    bed_file = tmp_path / "empty.bed"
    bed_file.write_text("chr3\t1000\t2000\n")
    bed = hl.import_bed(str(bed_file), reference_genome="GRCh38")
    result = filtering.filter_by_bed(mock_hail_table, bed)
    assert result.count() == 0
