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


def mutation_spectra_stats(df: pd.DataFrame, iqr_multiplier: float) -> (pd.DataFrame, pd.DataFrame):
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

    # Identify outliers
    is_outlier = df_melted.apply(
        lambda row: (row["Proportion"] < stats.at[row["mutation_type"], "lower"])
        | (row["Proportion"] > stats.at[row["mutation_type"], "upper"]),
        axis=1,
    )
    outliers = df_melted[is_outlier]

    return stats, outliers


def plot_mutation_spectra_boxplots(stats: pd.DataFrame, outliers: pd.DataFrame, width=800, height=600, **kwargs):
    """
    Plot mutation spectra by calculated table
    """
    mutation_types = stats.index.tolist()
    # Create the source for plotting
    source_data = {
        "mutation_type": mutation_types,
        "median": stats["50%"],
        "q25": stats["25%"],
        "q75": stats["75%"],
        "upper": stats["upper"],
        "lower": stats["lower"],
    }
    source = bokeh.models.ColumnDataSource(source_data)

    # Create the figure
    p = bokeh.plotting.figure(
        width=width, height=height, x_range=mutation_types, title="Mutation Spectrum", toolbar_location="above"
    )

    # quantile boxes
    boxplot_glyphs = p.vbar(
        x="mutation_type", width=0.7, bottom="q25", top="q75", source=source, line_color="black", fill_color="skyblue"
    )
    # Median line
    p.vbar(
        x="mutation_type",
        width=0.7,
        bottom="median",
        top="median",
        source=source,
        line_color="black",
        fill_color="skyblue",
    )

    boxplot_hover = bokeh.models.HoverTool(
        renderers=[boxplot_glyphs],
        tooltips=[
            ("Mutation Type", "@mutation_type"),
            ("Median", "@median{0.000}"),
            ("Q25", "@q25{0.000}"),
            ("Q75", "@q75{0.000}"),
        ],
    )
    p.add_tools(boxplot_hover)

    # outlier range
    whisker = bokeh.models.Whisker(base="mutation_type", upper="upper", lower="lower", source=source)
    whisker.upper_head.size = whisker.lower_head.size = 20
    p.add_layout(whisker)

    # Plot outliers
    outlier_glyphs = p.scatter("mutation_type", "Proportion", source=outliers, size=6, color="black", alpha=0.3)
    outlier_hover = bokeh.models.HoverTool(
        renderers=[outlier_glyphs], tooltips=[("Sample", "@samples"), ("Proportion", "@Proportion{0.0000}")]
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


def plot_mutation_spectra_barplot(stats: pd.DataFrame, width=800, height=600, **kwargs):
    """
    Plot mutation spectra from mutation spectra stats as a bar plot
    """
    # Define color mapping for symmetric pairs in mutation spectra
    MUTATION_SPECTRA_COLOR_MAPPING: dict[str, str] = {
        # A<->T base pairs
        "A>C": "#1f77b4",  # blue
        "T>G": "#1f77b4",  # blue (symmetric with A>C)
        "A>G": "#2ca02c",  # green
        "T>C": "#2ca02c",  # green (symmetric with A>G)
        "A>T": "#d62728",  # red
        "T>A": "#d62728",  # red (symmetric with A>T)
        # G<->C base pairs
        "G>A": "#ff7f0e",  # orange
        "C>T": "#ff7f0e",  # orange (symmetric with G>A)
        "G>C": "#9467bd",  # purple
        "C>G": "#9467bd",  # purple (symmetric with G>C)
        "G>T": "#8c564b",  # brown
        "C>A": "#8c564b",  # brown (symmetric with G>T)
    }
    # Color for any non-mapped SNVs. Should never appear in correct data.
    MUT_SPECTRA_NON_MAPPED_COLOR: str = "black"

    mutation_types = stats.index.tolist()
    # Create the source for plotting
    source_data = {
        "mutation_type": mutation_types,
        "median": stats["50%"],
        "std": stats["std"],
        "q25": stats["25%"],
        "q75": stats["75%"],
        "color": [
            MUTATION_SPECTRA_COLOR_MAPPING.get(mut, MUT_SPECTRA_NON_MAPPED_COLOR) for mut in mutation_types
        ],  # default to pink if type not found
    }
    source = bokeh.models.ColumnDataSource(source_data)

    # Create the figure
    p = bokeh.plotting.figure(
        width=width, height=height, x_range=mutation_types, title="Mutation Spectrum", toolbar_location="above"
    )

    # Plot median values as bars
    p.vbar(x="mutation_type", top="median", width=0.5, source=source, fill_color="color", line_color="black")

    # Add error whiskers for quartiles
    error_whisker = bokeh.models.Whisker(
        base="mutation_type", upper="std", lower="std", source=source, line_color="black"
    )
    error_whisker.upper_head.size = error_whisker.lower_head.size = 10
    p.add_layout(error_whisker)

    # Customize the plot
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = "gray"
    p.ygrid.grid_line_alpha = 0.1
    p.xaxis.axis_label = "Mutation Type"
    p.yaxis.axis_label = "Proportion"
    p.xaxis.major_label_orientation = math.pi / 4

    # Add hover tool
    hover = bokeh.models.HoverTool(
        tooltips=[
            ("Mutation Type", "@mutation_type"),
            ("Median", "@median{0.000}"),
            ("Standard Deviation", "@std{0.000}"),
            ("Q25", "@q25{0.000}"),
            ("Q75", "@q75{0.000}"),
        ]
    )
    p.add_tools(hover)

    return p


def plot_mutation_spectra_combined(mut_spectra, iqr_multiplier: float = 1.5, **kwargs):
    stats, outliers = mutation_spectra_stats(mut_spectra, iqr_multiplier)
    # Make and save the plot

    p_boxplot = plot_mutation_spectra_boxplots(stats, outliers, **kwargs)
    p_barplot = plot_mutation_spectra_barplot(stats, **kwargs)
    # Create two panels, one for each plot
    tab_boxplot = bokeh.models.TabPanel(child=p_boxplot, title="Boxplot")
    tab_barplot = bokeh.models.TabPanel(child=p_barplot, title="Barplot")
    # Combine the panels into mut_spectra_plot
    mut_spectra_plot = bokeh.models.Tabs(tabs=[tab_boxplot, tab_barplot])
    return mut_spectra_plot
