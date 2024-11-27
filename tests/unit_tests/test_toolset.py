import pathlib
import hail as hl  # type: ignore
from pytest import mark as m
from wes_qc import hail_utils, filtering, visualize


@m.context("When the temp folder does not have file:// prefix")
@m.it("should return None and do nothing")
def test_clear_temp_folder_1(tmp_path: pathlib.Path) -> None:
    tempdir = tmp_path / "temp"
    tempfile = tempdir / "test.txt"
    tempdir.mkdir()
    tempfile.touch()
    hail_utils.clear_temp_folder(str(tempdir))
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


@m.context("Plots results of the population PCA")
@m.it("Should run without exceptions")
def test_plot_pop_pca(pca_scores_table, tmp_path):
    plot_path = tmp_path / "test_pop_pca_plot.html"
    visualize.plot_pop_pca(pca_scores_table, str(plot_path), n_pca=3, pop="known_pop")
    assert plot_path.exists()
