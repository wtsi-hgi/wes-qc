import pathlib
from pytest import mark as m
from wes_qc import hail_utils, visualize


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
@m.it("Tests that hail inits without errors and can be stopped without errors")
def test_hail_init(tmp_path):
    hail_utils.init_hl(str(tmp_path))


@m.context("Plots results of the population PCA")
@m.it("Should run without exceptions")
def test_plot_pop_pca(pca_scores_table, tmp_path):
    plot_path = tmp_path / "test_pop_pca_plot.html"
    visualize.plot_pop_pca(pca_scores_table, str(plot_path), n_pca=3, pop="known_pop")
    assert plot_path.exists()
    plot_path = tmp_path / "test_pop_pca_plot_no_pop_color.html"
    visualize.plot_pop_pca(pca_scores_table, str(plot_path), n_pca=3, pop=None)
    assert plot_path.exists()
