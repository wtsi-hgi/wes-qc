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

    # Make and save the plot
    p_box_whiskers = visualize.plot_mutation_spectra(mut_spectra, **config["step1"]["plot_mutation_spectra_preqc"])
    p_hist = visualize.plot_simplified_mutation_spectra(
        df=mut_spectra, **config["step1"]["plot_mutation_spectra_preqc"]
    )

    from bokeh.models import Tabs, TabPanel

    # Assuming p_box_whiskers and p_hist are your Bokeh plot objects

    # Create two panels, one for each plot
    tab1 = TabPanel(child=p_box_whiskers, title="Boxplot")
    tab2 = TabPanel(child=p_hist, title="Histogram")

    # Combine the panels into tabs
    tabs = Tabs(tabs=[tab1, tab2])

    # Save the combined plot
    bokeh.io.output_file(mut_spectra_plot_path)
    bokeh.io.save(tabs)

    # bokeh.io.output_file(mut_spectra_plot_path)
    # bokeh.io.save(p_box_whiskers)


if __name__ == "__main__":
    main()
