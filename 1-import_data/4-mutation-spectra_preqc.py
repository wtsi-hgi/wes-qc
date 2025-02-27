import bokeh.io
import hail as hl

from utils.utils import parse_config, path_spark
from wes_qc import hail_utils, visualize


def main():
    # set up
    config = parse_config()
    tmp_dir = config["general"]["tmp_dir"]

    # = STEP DEPENDENCIES = #
    mtpath = config["step1"]["gatk_mt_outfile"]

    # = STEP OUTPUTS = #
    mut_spectra_path = config["step1"]["plot_mutation_spectra_preqc"]["mut_spectra_path"]
    mut_spectra_plot_path = config["step1"]["plot_mutation_spectra_preqc"]["mut_spectra_plot_path"]

    # = STEP LOGIC = #
    hail_utils.init_hl(tmp_dir)

    # Calculate mutation spectra
    mt = hl.read_matrix_table(path_spark(mtpath))
    mut_spectra = visualize.calculate_mutation_spectra(mt)
    mut_spectra.to_csv(mut_spectra_path, sep="\t")

    # Save the combined plot
    mut_spectra_plot = visualize.plot_mutation_spectra_combined(
        mut_spectra, **config["step1"]["plot_mutation_spectra_preqc"]
    )
    bokeh.io.output_file(mut_spectra_plot_path)
    bokeh.io.save(mut_spectra_plot)


if __name__ == "__main__":
    main()
