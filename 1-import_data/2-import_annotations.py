import hail as hl  # type: ignore
import hail.plot
import bokeh

from utils.utils import parse_config
from utils.config import path_spark
from wes_qc import hail_utils

verifybamid_types = {
    "#SEQ_ID": "str",
    "RG": "str",
    "CHIP_ID": "str",
    "#SNPS": "int",
    "#READS": "int",
    "AVG_DP": "float",
    "FREEMIX": "float",
    "FREELK1": "float",
    "FREELK0": "float",
    "FREE_RH": "str",
    "FREE_RA": "str",
    "CHIPMIX": "str",
    "CHIPLK1": "str",
    "CHIPLK0": "str",
    "CHIP_RH": "str",
    "CHIP_RA": "str",
    "DPREF": "str",
    "RDPHET": "str",
    "RDPALT": "str",
}


def validate_verifybamid(
    mt: hl.MatrixTable,
    verifybamid: hl.Table,
    verifybamid_plot: str,
    samples_failing_freemix_tsv: str,
    freemix_treshold: float = 0.05,
    **kwargs,
) -> hl.MatrixTable:
    verifybamid = verifybamid.key_by("#SEQ_ID")
    mt = mt.annotate_cols(freemix=verifybamid[mt.s].FREEMIX)

    samples = mt.cols()
    samples = samples.order_by(hl.asc(samples.freemix))
    samples = samples.add_index(name="s_index")

    # Exporting the list of samples failing Freemix
    samples_failing_freemix = samples.filter(samples.freemix > freemix_treshold)
    n_samples_failed_freemix = samples_failing_freemix.count()
    if n_samples_failed_freemix > 0:
        samples_failing_freemix.export(path_spark(samples_failing_freemix_tsv))
        print(f"=== {n_samples_failed_freemix} samples failed freemix")
        print(f"=== Failing samples exported to {samples_failing_freemix_tsv}")

    samples_without_freemix = samples.filter(~hl.is_defined(samples.freemix))
    n_samples_without_freemix = samples_without_freemix.count()
    if n_samples_without_freemix > 0:
        print(f"=== Detected {n_samples_without_freemix} samples without freemix: ", end="")
        print(" ".join(samples_without_freemix.s.collect()))

    # Plottign freemix score
    p = hail.plot.scatter(
        samples.s_index,
        samples.freemix,
        xlabel="n",
        ylabel="Freemix score",
        title="Freemix score",
        width=800,
        height=800,
    )
    hline = bokeh.models.Span(location=freemix_treshold, dimension="width", line_color="red", line_width=3)
    # Add the Span annotation to the plot
    p.renderers.extend([hline])

    label = bokeh.models.Label(
        x=20,
        y=700,
        text=f"{n_samples_failed_freemix} samples failed freemix",
        x_units="screen",
        y_units="screen",
        text_color="black",
        text_font_size="14pt",
    )
    p.add_layout(label)

    bokeh.plotting.output_file(verifybamid_plot)
    bokeh.plotting.save(p)

    return mt


def main() -> None:
    # = STEP SETUP = #
    config = parse_config()

    tmp_dir = config["general"]["tmp_dir"]

    # = STEP PARAMETERS = #

    # = STEP DEPENDENCIES = #
    mtpath = config["step1"]["gatk_mt_outfile"]
    verifybamid_selfsm = config["step1"]["validate_verifybamid"]["verifybamid_selfsm"]

    # = STEP OUTPUTS = #
    mtoutpath = config["step1"]["validate_verifybamid"]["mt_with_freemix"]

    # = STEP LOGIC = #
    _ = hail_utils.init_hl(tmp_dir)

    mt = hl.read_matrix_table(path_spark(mtpath))
    verifybamid = hl.import_table(path_spark(verifybamid_selfsm), types=verifybamid_types)
    mt = validate_verifybamid(mt, verifybamid, **config["step1"]["validate_verifybamid"])
    mt.write(path_spark(mtoutpath), overwrite=True)


if __name__ == "__main__":
    main()
