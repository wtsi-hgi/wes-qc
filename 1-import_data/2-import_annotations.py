import hail as hl  # type: ignore
import hail.plot
import bokeh

from utils.utils import parse_config
from wes_qc.hail_utils import path_spark
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
    else:
        print("=== OK: All samples passed Freemix validation")

    samples_without_freemix = samples.filter(~hl.is_defined(samples.freemix))
    n_samples_without_freemix = samples_without_freemix.count()
    if n_samples_without_freemix > 0:
        print(f"=== WARNING: Detected {n_samples_without_freemix} samples without freemix: ", end="")
        print(" ".join(samples_without_freemix.s.collect()))
    else:
        print("=== OK: All samples are annotated with FreeMix scores")

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


def annotate_self_reported_sex(mt: hl.MatrixTable, sex_metadata_file: str, **kwargs) -> hl.MatrixTable:
    """
    Annotates samples in the matrix-table with the self-reported sex
    """
    metadata_ht = hl.import_table(path_spark(sex_metadata_file), delimiter="\t").key_by("sample_id")
    metadata_ht = metadata_ht.transmute(self_reported_sex=metadata_ht.self_reported_sex.lower())
    mt_sex_annotated = mt.annotate_cols(self_reported_sex=metadata_ht[mt.s].self_reported_sex)

    # Checking for the samples without self-reported sex
    samples = mt_sex_annotated.cols()
    samples_without_self_reported_sex = samples.filter(~hl.is_defined(samples.self_reported_sex))
    n_samples_no_self_reported_sex = samples_without_self_reported_sex.count()
    if n_samples_no_self_reported_sex > 0:
        print(f"=== WARNING: Detected {n_samples_no_self_reported_sex} samples without self-reported sex: ", end="")
        print(" ".join(samples_without_self_reported_sex.s.collect()))
    else:
        print("=== OK: All samples are annotated with self-reported sex")

    return mt_sex_annotated


def main() -> None:
    # = STEP SETUP = #
    config = parse_config()

    tmp_dir = config["general"]["tmp_dir"]

    # = STEP PARAMETERS = #

    # = STEP DEPENDENCIES = #
    mtpath = config["step1"]["gatk_mt_outfile"]
    verifybamid_selfsm = config["step1"]["validate_verifybamid"]["verifybamid_selfsm"]
    # TODO: change after merging previous config changes from main
    sex_metadata_file = config["step1"]["sex_metadata_file"]

    # = STEP OUTPUTS = #
    mtoutpath = config["step1"]["mt_metadata_annotated"]

    # = STEP LOGIC = #
    _ = hail_utils.init_hl(tmp_dir)

    mt = hl.read_matrix_table(path_spark(mtpath))

    if verifybamid_selfsm is not None:
        print("=== Running verifyBamID validation ")
        verifybamid = hl.import_table(path_spark(verifybamid_selfsm), types=verifybamid_types)
        mt = validate_verifybamid(mt, verifybamid, **config["step1"]["validate_verifybamid"])
    else:
        print("=== Skipping verifyBamID validation")

    if sex_metadata_file is not None:
        print("=== Annotating self-reported sex ")
        mt = annotate_self_reported_sex(mt, path_spark(sex_metadata_file))
    else:
        print("=== Skipping self-reported sex annotation")

    mt.write(path_spark(mtoutpath), overwrite=True)


if __name__ == "__main__":
    main()
