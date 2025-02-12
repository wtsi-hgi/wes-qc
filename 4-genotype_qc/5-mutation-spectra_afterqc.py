import bokeh.io
import hail as hl

from utils.utils import parse_config, path_spark
from wes_qc import hail_utils, visualize


def main():
    # set up
    config = parse_config()
    tmp_dir = config["general"]["tmp_dir"]

    # = STEP DEPENDENCIES = #
    mtpath = config["step4"]["export_vcfs_a"]["mtfile"]

    # = STEP OUTPUTS = #
    mut_spectra_path = config["step4"]["plot_mutation_spectra_afterqc"]["mut_spectra_path"]
    mut_spectra_plot_path = config["step4"]["plot_mutation_spectra_afterqc"]["mut_spectra_plot_path"]

    # = STEP LOGIC = #
    hail_utils.init_hl(tmp_dir)

    mt = hl.read_matrix_table(path_spark(mtpath))
    mut_spectra = visualize.calculate_mutation_spectra(mt)
    mut_spectra.to_csv(mut_spectra_path, sep="\t")

    # Save the plot
    p = visualize.plot_mutation_spectra(mut_spectra, **config["step4"]["plot_mutation_spectra_afterqc"])
    bokeh.io.output_file(mut_spectra_plot_path)
    bokeh.io.save(p)


if __name__ == "__main__":
    main()
