# perform hail sample QC stratified by superpopulation and identify outliers
import os.path

import hail as hl
import pandas as pd

from wes_qc import hail_utils
from utils.utils import parse_config, path_spark
from gnomad.sample_qc.filtering import compute_stratified_metrics_filter
import bokeh.plotting as bkplot
import bokeh.layouts as bklayouts
from bokeh.models import Div, Span, Range1d, Label, Title
import numpy as np
from collections import defaultdict


# TODO: rename to annotate_with_pop
def annotate_mt(raw_mt_file: str, pop_ht_file: str):
    """
    Annotate mt with superpopulation and sequencing runid
    :param str raw_mt_file: raw mt file
    :param str pop_ht_file: file with population annotations
    :param str runid_file: metadata file to annotate run ids TODO
    :param str annotated_mt_file: annotated mt file
    :param str pop_pandas_file: tsv from pandas df - EGA and pop
    :param dict config: A config object. No effect.

    ### Config fields
    None

    ### Indirect config fields
    step1.gatk_mt_outfile : input path : used in main
    step2.predict_populations.pop_ht_file : input path : used in main
    step2.annotate_with_pop.annotated_mt_file : output path : used in main
    """
    mt = hl.read_matrix_table(path_spark(raw_mt_file))
    pop_ht = hl.read_table(path_spark(pop_ht_file))
    mt = mt.annotate_cols(assigned_pop=pop_ht[mt.s].pop)
    return mt


def stratified_sample_qc(
    mt: hl.MatrixTable,
    min_depth: int,
    min_genotype_quality: float,
    min_vaf: float,
    compute_stratified_metrics_filter_args: dict,
    **kwargs,
):
    """
    Run sample QC and stratify by population
    param str annotated_mt_file: population and run id annotated MT file
    param str mt_qc_outfile: sample QC MT file
    param str ht_qccols_outfile: sample QC columns HT file
    param str qc_filter_file: output file for stratified sample QC HT
    param dict config: config object

    TODO: note about `if min_depth > 0 or min_genotype_quality > 0 or min_vaf > 0`

    ### Config fields
    step2.stratified_sample_qc.min_depth : float : TODO
    step2.stratified_sample_qc.min_genotype_quality : float : TODO
    step2.stratified_sample_qc.min_vaf : float : TODO
    step2.stratified_sample_qc.output_text_file : output path : TODO
    step2.stratified_sample_qc.output_globals_json_file : output path : TODO

    ### Indirect config fields
    step2.annotate_with_pop.annotated_mt_file : input path : used in main
    step2.stratified_sample_qc.mt_qc_outfile : output path : used in main
    step2.stratified_sample_qc.ht_qc_cols_outfile : output path : used in main
    step2.stratified_sample_qc.qc_filter_file : output path : used in main
    """
    # filter to autosomes only
    mt = mt.filter_rows(mt.locus.in_autosome())

    # filter MT by depth/gq/vaf
    if min_depth > 0 or min_genotype_quality > 0 or min_vaf > 0:
        vaf = mt.AD[1] / hl.sum(mt.AD)
        print(
            f"Filtering input MT by depth: DP={min_depth}, genotype quality: GQ={min_genotype_quality}, VAF: VAF={min_vaf}"
        )
        filter_condition = (
            (mt.GT.is_het() & (vaf > min_vaf) & (mt.DP > min_depth) & (mt.GQ > min_genotype_quality))
            | (mt.GT.is_hom_ref() & (mt.DP > min_depth) & (mt.GQ > min_genotype_quality))
            | (mt.GT.is_hom_var() & (mt.DP > min_depth) & (mt.GQ > min_genotype_quality))
        )
        fraction_filtered = mt.aggregate_entries(hl.agg.fraction(~filter_condition))
        print(f"Filtering {fraction_filtered * 100:.2f}% entries out of downstream analysis.")
        mt = mt.filter_entries(filter_condition)

    # run sample QC
    mt_with_sampleqc = hl.sample_qc(mt, name="sample_qc")
    # annotate with heterozygosity rate
    mt_with_sampleqc = mt_with_sampleqc.annotate_cols(
        sample_qc=mt_with_sampleqc.sample_qc.annotate(
            heterozygosity_rate=mt_with_sampleqc.sample_qc.n_het / mt_with_sampleqc.sample_qc.n_called
        )
    )

    # stratify by pop
    pop_ht = mt_with_sampleqc.cols()
    # TODO: shall we move qc_metrics to the config?
    qc_metrics = [
        "heterozygosity_rate",
        "n_snp",
        "r_ti_tv",
        "n_transition",
        "n_transversion",
        "r_insertion_deletion",
        "n_insertion",
        "n_deletion",
        "r_het_hom_var",
    ]

    # Using gnomAD function to calculate stratified metrics
    pop_filter_ht = compute_stratified_metrics_filter(
        pop_ht,
        qc_metrics={metric: pop_ht.sample_qc[metric] for metric in qc_metrics},
        strata={"qc_pop": pop_ht.assigned_pop},
        **compute_stratified_metrics_filter_args,
    )

    globals = hl.eval(pop_filter_ht.globals.qc_metrics_stats)

    pop_ht = pop_ht.annotate_globals(qc_metrics_stats=globals)
    pop_ht = pop_ht.annotate(**pop_filter_ht[pop_ht.key]).persist()
    checkpoint = pop_ht.aggregate(hl.agg.count_where(hl.len(pop_ht.qc_metrics_filters) == 0))
    print(f"=== Samples passing pop filtering: {checkpoint}")

    # return the pop_ht for plotting in the next step
    return mt_with_sampleqc, pop_ht


def add_caption_to_plot(p, metric, pop):
    caption_title = Title(text=f"{metric} for {pop}", text_font_size="14pt", align="center")
    p.add_layout(caption_title, "above")
    p.yaxis.axis_label = "Count"
    p.xaxis.axis_label = metric
    return p


def plot_sampleqc_metric(
    df: pd.Series,
    lower_threshold=None,
    upper_threshold=None,
    n_bins=150,
    color="blue",
    plot_width: int = 800,
    plot_height: int = 600,
):
    """
    Plots a histogram of the given metric values for a specific population, with options to add
    thresholds, annotations, captions, and configure visual aspects of the plot.

    Parameters:
        df: Input data series for the metric to be visualized.
        Should be already filtered to contain only single pop and a single metric
        metric: Metric name to be displayed on the x-axis and caption of the plot.
        pop: Name of the population to be included in the plot's caption.
        lower_threshold: A numeric value to mark the lower boundary on the plot.
           Default is None.
        upper_threshold: A numeric value to mark the upper boundary on the plot.
           Default is None.
        n_bins: Number of bins to use for the histogram.
        Default is 150.
        color: Color used for the histogram and plot elements.
        Default is 'blue'.
        plot_width: Width of the plot.
        Default is 800.
        plot_height: Height of the plot.
        Default is 600.

    Returns:
        bokeh.plotting.figure: A Bokeh plot object that visualizes the specified metric.
    """

    p = bkplot.figure(width=plot_width, height=plot_height)
    hist, edges = np.histogram(df[df.notna()], bins=n_bins)
    interval = edges[1] - edges[0]
    padding = max(1, n_bins // 10) * interval
    p.x_range.start = edges[0] - padding
    p.x_range.end = edges[-1] + padding

    p.quad(
        top=hist,
        bottom=0,
        left=edges[:-1],
        right=edges[1:],
        fill_color=color,
        line_color=color,
    )

    if lower_threshold is not None:
        n_below = len(df[df < lower_threshold])
        hline_lower = Span(location=lower_threshold, dimension="height", line_color="red", line_width=2)
        ann_below = Label(
            x=lower_threshold,
            y=hist.max(),
            text=f"{n_below}",
            text_font_size="10pt",
            background_fill_color="white",
            background_fill_alpha=0.5,
            text_align="right",
            x_offset=-5,
        )
        p.add_layout(hline_lower)
        p.add_layout(ann_below)

    if upper_threshold is not None:
        n_above = len(df[df > upper_threshold])
        hline_upper = Span(location=upper_threshold, dimension="height", line_color="red", line_width=2)

        ann_above = Label(
            x=upper_threshold,
            y=hist.max(),
            text=f"{n_above}",
            text_font_size="10pt",
            background_fill_color="white",
            background_fill_alpha=0.5,
            text_align="left",
            x_offset=5,
        )
        p.add_layout(hline_upper)
        p.add_layout(ann_above)
    return p


def plot_sample_qc_metrics(
    pop_ht: hl.Table, plot_outdir: str, plot_width: int = 400, plot_height: int = 400, n_bins: int = 100, **kwargs
):
    def extract_lower_upper(metadata):
        output = {}
        stat = metadata["qc_metrics_stats"]
        for k, v in stat.items():
            for k2, v2 in v.items():
                output[(k2, k[0])] = (v2["lower"], v2["upper"])
        return output

    os.makedirs(plot_outdir, exist_ok=True)

    pd_ht = pop_ht.to_pandas()
    stats = hl.eval(pop_ht.globals)

    limits = extract_lower_upper(stats)

    metrics = list(set([met[0] for met in limits.keys()]))

    pops = set()
    plots = []
    plots_by_metric = defaultdict(list)
    all_pops = ["EUR", "EAS", "AFR", "AMR", "SAS", "oth"]
    colors = ["#F0E442", "#D55E00", "#0072B2", "#009E73", "#E69F00", "#56B4E9"]
    pop_colors = dict(zip(all_pops, colors))

    for pop, subdf in pd_ht.groupby("assigned_pop"):
        pops.add(pop)
        for metric in metrics:
            data = subdf[f"sample_qc.{metric}"]
            lower, upper = limits[(metric, pop)]
            p = plot_sampleqc_metric(
                data,
                lower_threshold=lower,
                upper_threshold=upper,
                n_bins=n_bins,
                color=pop_colors[pop],
                plot_width=plot_width,
                plot_height=plot_height,
            )

            plots.append(p)
            plots_by_metric[metric].append(p)

            # Saving individual plot.
            # Need to create another plot due to Bokeh restrictions
            # Python deepcopy doesn't work
            p = plot_sampleqc_metric(
                data,
                lower_threshold=lower,
                upper_threshold=upper,
                n_bins=n_bins,
                color=pop_colors[pop],
                plot_width=plot_width,
                plot_height=plot_height,
            )
            plot_with_caption = add_caption_to_plot(p, metric, pop)
            plot_name = f"SampleQC_hist_{metric}_{pop}.html"
            bkplot.output_file(os.path.join(plot_outdir, plot_name))
            bkplot.save(plot_with_caption)

        plots.append(
            Div(
                text=f"<b>{pop}</b>",
                width=10,
                height=plot_height,
                styles={"writing-mode": "vertical-rl", "text-align": "center"},
            )
        )

    # select the minimum and maximum x range for each metric and share it across columns
    for k, p in plots_by_metric.items():
        x_lower = min([pl.x_range.start for pl in p])
        x_upper = max([pl.x_range.end for pl in p])
        for pl in p:
            pl.x_range = Range1d(start=x_lower, end=x_upper)

    col_titles = [
        Div(text=f"<b>{metric}</b>", width=plot_width, height=15, styles={"text-align": "center"}) for metric in metrics
    ] + [None]

    grid = bklayouts.gridplot(col_titles + plots, ncols=len(metrics) + 1)
    plot_outfile = os.path.join(plot_outdir, "sample_qc_all_metrics_by_pop.html")
    bkplot.output_file(plot_outfile)
    bkplot.save(grid)


def main():
    # = STEP SETUP = #
    config = parse_config()
    tmp_dir = config["general"]["tmp_dir"]

    # = STEP PARAMETERS = #

    # = STEP DEPENDENCIES = #
    raw_mt_file = config["step1"]["gatk_mt_outfile"]
    pop_ht_file = config["step2"]["predict_pops"]["pop_ht_outfile"]

    # = STEP OUTPUTS = #
    annotated_mt_file = config["step2"]["annotate_with_pop"]["annotated_mt_file"]
    mt_qc_outfile = config["step2"]["stratified_sample_qc"]["mt_qc_outfile"]
    ht_qc_cols_outfile = config["step2"]["stratified_sample_qc"]["ht_qc_cols_outfile"]
    qc_filter_file = config["step2"]["stratified_sample_qc"]["qc_filter_file"]
    output_text_file = config["step2"]["stratified_sample_qc"]["output_text_file"]
    output_globals_json = config["step2"]["stratified_sample_qc"]["output_globals_json_file"]

    # = STEP LOGIC = #
    _ = hail_utils.init_hl(tmp_dir)

    # annotate mt with runid and pop
    mt = annotate_mt(raw_mt_file, pop_ht_file)
    mt.write(path_spark(annotated_mt_file), overwrite=True)

    # run sample QC and stratify by population
    mt_with_sampleqc, pop_ht = stratified_sample_qc(mt, **config["step2"]["stratified_sample_qc"])
    mt_with_sampleqc.write(path_spark(mt_qc_outfile), overwrite=True)
    mt_with_sampleqc.cols().write(path_spark(ht_qc_cols_outfile), overwrite=True)
    pop_ht.write(path_spark(qc_filter_file), overwrite=True)
    pop_ht.export(path_spark(output_text_file), delimiter="\t")
    pop_ht.globals.export(path_spark(output_globals_json))

    # plot population metrics
    print("=== Plotting population metrics ===")
    plot_sample_qc_metrics(pop_ht, **config["step2"]["plot_sample_qc_metrics"])


if __name__ == "__main__":
    main()
