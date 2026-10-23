import hail as hl  # type: ignore
import hail.plot
import bokeh

from wxs_qc.config import get_config
from wxs_qc.hail_utils import path_spark
from wxs_qc import hail_utils, constants
from wxs_qc.input_validation import validate_self_reported_sex

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
    remove_freemix_outliers: bool = False,
    width=constants.width,
    height=constants.height,
    text_size=constants.plots_text_size,
    **kwargs,
) -> hl.MatrixTable:
    """
    Validates and visualizes Freemix scores from verifyBamID against a Hail MatrixTable.
    This function annotates the MatrixTable with Freemix scores, identifies samples failing the Freemix threshold,
    exports them if applicable, and generates a scatter plot to visualize the Freemix scores.
    """
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
        width=width,
        height=height,
    )
    hline = bokeh.models.Span(location=freemix_treshold, dimension="width", line_color="red", line_width=3)
    # Add the Span annotation to the plot
    p.renderers.extend([hline])

    label = bokeh.models.Label(
        x=20,
        y=height - 200,
        text=f"{n_samples_failed_freemix} samples failed freemix",
        x_units="screen",
        y_units="screen",
        text_color="black",
        text_font_size=text_size,
    )
    p.add_layout(label)

    p.axis.axis_label_text_font_size = text_size
    p.axis.major_label_text_font_size = text_size
    p.title.text_font_size = text_size

    bokeh.plotting.output_file(verifybamid_plot)
    bokeh.plotting.save(p)

    if remove_freemix_outliers:
        print("=== Removing freemix outliers")
        # Removing Freemix outliers
        mt = mt.filter_cols((mt.freemix <= freemix_treshold) | hl.is_missing(mt.freemix))
    return mt


def annotate_self_reported_sex(mt: hl.MatrixTable, sex_metadata_file: str, **kwargs) -> hl.MatrixTable:
    """
    Annotates samples in the matrix-table with the self-reported sex
    """
    metadata_ht = hl.import_table(path_spark(sex_metadata_file), delimiter="\t").key_by("sample_id")
    metadata_ht = metadata_ht.transmute(self_reported_sex=metadata_ht.self_reported_sex.lower())
    validate_self_reported_sex(metadata_ht)
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


def annotate_batch(mt: hl.MatrixTable, batch_metadata_file: str, **kwargs) -> hl.MatrixTable:
    """
    Annotates samples in the matrix-table with the batch information
    """
    metadata_ht = hl.import_table(path_spark(batch_metadata_file), delimiter="\t").key_by("sample_id")
    mt_batch_annotated = mt.annotate_cols(batch=metadata_ht[mt.s].batch)

    # Checking for the samples without self-reported sex
    samples = mt_batch_annotated.cols()
    samples_without_batch = samples.filter(~hl.is_defined(samples.batch))
    n_samples_no_batch = samples_without_batch.count()
    if n_samples_no_batch > 0:
        print(f"=== WARNING: Detected {n_samples_no_batch} samples without batch information: ", end="")
        print(" ".join(samples_without_batch.s.collect()))
    else:
        print("=== OK: All samples are annotated with batch information")
    return mt_batch_annotated


def main() -> None:
    # = STEP SETUP = #
    config = get_config()

    tmp_dir = config["general"]["tmp_dir"]

    # = STEP PARAMETERS = #

    # = STEP DEPENDENCIES = #
    mtpath = config["stage1"]["gatk_mt_outfile"]
    verifybamid_selfsm = config["stage1"]["validate_verifybamid"]["verifybamid_selfsm"]
    # TODO: change after merging previous config changes from main
    sex_metadata_file = config["stage1"]["sex_metadata_file"]
    batch_metadata_file = config["stage1"]["batch_metadata_file"]
    # = STEP OUTPUTS = #
    mtoutpath = config["stage1"]["mt_metadata_annotated"]

    # = STEP LOGIC = #
    _ = hail_utils.init_hl(tmp_dir)

    mt = hl.read_matrix_table(path_spark(mtpath))

    if verifybamid_selfsm is not None:
        print("=== Running verifyBamID validation ")
        verifybamid = hl.import_table(path_spark(verifybamid_selfsm), types=verifybamid_types)
        mt = validate_verifybamid(mt, verifybamid, **config["stage1"]["validate_verifybamid"])

    else:
        print("=== Skipping verifyBamID validation")

    if sex_metadata_file is not None:
        print("=== Annotating self-reported sex ")
        mt = annotate_self_reported_sex(mt, path_spark(sex_metadata_file))
    else:
        print("=== Skipping self-reported sex annotation")

    if batch_metadata_file is not None:
        print("=== Annotating batch ")
        mt = annotate_batch(mt, path_spark(batch_metadata_file))
    else:
        print("=== Skipping batch annotation")

    mt.write(path_spark(mtoutpath), overwrite=True)


if __name__ == "__main__":
    main()
