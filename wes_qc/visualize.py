from typing import Optional

import hail as hl
import bokeh
import math

import pandas as pd


def plot_pop_pca(pca_scores: hl.Table, plot_file: str, n_pca: int = 3, pop: Optional[str] = None) -> None:
    """
    Plots grid of n_pca PCA components
    """
    # Making fixed color mapping for superpopulations
    pops = ["EUR", "EAS", "AFR", "AMR", "SAS", "oth"]
    colors = ["#F0E442", "#D55E00", "#0072B2", "#009E73", "#E69F00", "#56B4E9"]
    pop_colors_mapper = bokeh.models.CategoricalColorMapper(factors=pops, palette=colors)

    colors_map = pop_colors_mapper if pop is not None else None

    layout = [[None] * n_pca for i in range(n_pca)]
    label = pca_scores[pop] if pop is not None else None
    for i in range(n_pca):
        for j in range(i + 1, n_pca):
            p = hl.plot.scatter(
                pca_scores.scores[i],
                pca_scores.scores[j],
                xlabel=f"PC{i+1}",
                ylabel=f"PC{j+1}",
                label=label,
                colors=colors_map,
                title=f"PC{i+1} vs PC{j+1}",
            )
            layout[i][j] = p

    plots = bokeh.layouts.gridplot(layout)
    bokeh.plotting.output_file(plot_file)
    bokeh.plotting.save(plots)


def calculate_mutation_spectra(mt):
    """
    Calculates the mutation spectra from the matritable and returns it in the pivoted Pandas dataframe
    """
    # calculate mutation fraction
    mt = mt.filter_rows(mt.locus.in_autosome())

    mt = mt.filter_rows((hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1]))

    mt = mt.annotate_rows(mutation_type=hl.str(mt.alleles[0]) + ">" + hl.str(mt.alleles[1]))

    grouped_mutations = mt.group_rows_by(mt.mutation_type).aggregate(
        alt_call_n=hl.agg.filter((mt.GT.is_het()), hl.agg.count())
        + 2 * (hl.agg.filter((mt.GT.is_hom_var()), hl.agg.count()))
    )

    mutation_df = (
        grouped_mutations.entries()
        .select("alt_call_n")
        .to_pandas()
        .pivot(index="s", columns="mutation_type", values="alt_call_n")
    )
    fraction_mutation_df = mutation_df.div(mutation_df.sum(axis=1), axis=0)
    return fraction_mutation_df


def plot_mutation_spectra(df: pd.DataFrame, iqr_multiplier: float = 1.5, width=800, height=600, **kwargs):
    """
    Plot mutation spectra by calculated table
    """
    # Prepare the data
    df_copy = df.copy()
    df_copy["samples"] = df_copy.index
    df_melted = df_copy.melt(id_vars=["samples"], var_name="mutation_type", value_name="Proportion")
    df_melted["Proportion"] = df_melted["Proportion"].astype(
        "float64"
    )  # For some samples Hail returns NaN that Pandas do not recognize

    # Calculate statistics for box plot
    stats = df_melted.groupby("mutation_type").Proportion.describe()
    iqr = stats["75%"] - stats["25%"]
    stats["lower"] = stats["25%"] - iqr_multiplier * iqr
    stats["upper"] = stats["75%"] + iqr_multiplier * iqr

    mutation_types = stats.index.tolist()
    source = bokeh.models.ColumnDataSource(stats)

    # Create the figure
    p = bokeh.plotting.figure(
        width=width, height=height, x_range=mutation_types, title="Mutation Spectrum", toolbar_location="above"
    )

    # quantile boxes
    p.vbar(
        x="mutation_type", width=0.7, bottom="25%", top="50%", source=source, line_color="black", fill_color="skyblue"
    )
    p.vbar(
        x="mutation_type", width=0.7, bottom="50%", top="75%", source=source, line_color="black", fill_color="skyblue"
    )

    # outlier range
    whisker = bokeh.models.Whisker(base="mutation_type", upper="upper", lower="lower", source=source)
    whisker.upper_head.size = whisker.lower_head.size = 20
    p.add_layout(whisker)

    # Identify outliers
    is_outlier = df_melted.apply(
        lambda row: (row["Proportion"] < stats.at[row["mutation_type"], "lower"])
        | (row["Proportion"] > stats.at[row["mutation_type"], "upper"]),
        axis=1,
    )
    outliers = df_melted[is_outlier]
    # Plot outliers
    outlier_glyphs = p.scatter("mutation_type", "Proportion", source=outliers, size=6, color="black", alpha=0.3)

    outlier_hover = bokeh.models.HoverTool(
        renderers=[outlier_glyphs], tooltips=[("Sample", "@samples"), ("Proportion", "@Proportion")]
    )
    p.add_tools(outlier_hover)

    # Customize the plot
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = "gray"
    p.ygrid.grid_line_alpha = 0.1
    p.xaxis.axis_label = "Mutation Type"
    p.yaxis.axis_label = "Proportion"
    p.xaxis.major_label_orientation = math.pi / 4

    return p
