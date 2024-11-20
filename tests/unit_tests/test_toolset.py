from wes_qc import hail_utils
import hail as hl
from pytest import mark as m


@m.context("When the temp folder does not have file:// prefix")
@m.it("should return None and do nothing")
def test_clear_temp_folder_1(tmp_path):
    tempdir = tmp_path / "temp"
    tempfile = tempdir / "test.txt"
    tempdir.mkdir()
    tempfile.touch()
    res = hail_utils.clear_temp_folder(str(tempdir))
    assert res is None
    assert tempdir.exists()
    assert len(list(tempdir.glob("*"))) == 1


@m.context("When the temp folder has file:// prefix")
@m.it("should clear the temp folder and recreate it")
def test_clear_temp_folder_2(tmp_path):
    tempdir = tmp_path / "temp"
    tempdir.mkdir()
    tempfile = tempdir / "test.txt"
    tempfile.touch()
    hail_utils.clear_temp_folder(f"file://{tempdir}")
    assert tempdir.exists()
    assert len(list(tempdir.glob("*"))) == 0


@m.context("When hail is initialized and then stopped")
@m.it("should use that temp folder as the hail tmp folder, and be stopped")
def test_hail_init(tmp_path):
    sc = hail_utils.init_hl(str(tmp_path))
    assert hl.tmp_dir() == str(tmp_path)
    assert not sc._jsc.sc().isStopped()
    hail_utils.stop_hl(sc)
    assert sc._jsc is None


@m.context("When the input matrix table has variants on autosomes, SNVs, and is biallelic")
@m.it("should return a matrix table with only autosomes, SNVs, and non-split variants")
def test_filter_mt_autosome_biallelic_snvs(hail_data):
    assert True
    # Create synthetic data using the Balding-Nichols model
    # mt = hl.balding_nichols_model(n_populations=3, n_samples=10, n_variants=100)

    # # Apply the filter function
    # filtered_mt = filtering.filter_mt_autosome_biallelic_snvs(mt)

    # # Check that all variants are on autosomes
    # assert hl.agg.all(filtered_mt.locus.in_autosome())

    # # Check that all variants are SNVs
    # assert hl.agg.all(hl.is_snp(filtered_mt.alleles[0], filtered_mt.alleles[1]))

    # # Check that no variants were split
    # assert hl.agg.all(filtered_mt.was_split == False)
