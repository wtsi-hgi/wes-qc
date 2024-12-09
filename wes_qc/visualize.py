import os
from typing import Optional

import hail as hl
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import bokeh


def plot_pca_seaborn(pca_scores: hl.Table, plot_file: str, n_pca: int = 3, pop: Optional[str] = None) -> None:
    sns_data = pd.DataFrame.from_records(pca_scores.to_pandas().scores)[range(n_pca)]
    sns_data["s"] = pca_scores.s.collect()
    if pop is not None:
        sns_data[pop] = pca_scores[pop].collect()

    p = sns.pairplot(sns_data, vars=sns_data.columns[range(n_pca)], plot_kws={"s": 2}, hue=pop)
    p.savefig(plot_file)

    # p = seaborn.pairplot(sns_data1, vars=sns_data1.columns[range(5)], hue='SR_Ethnicity', plot_kws={'s': 2})


def plot_pca_bokeh(pca_scores: hl.Table, plot_file: str, n_pca: int = 3, pop: Optional[str] = None) -> None:
    """
    Plots grid of n_pca PCA components
    """
    # Making fixed color mapping for superpopulations
    pops = ["EUR", "EAS", "AFR", "AMR", "SAS", "oth"]
    colors = ["green", "goldenrod", "brown", "indigo", "red", "grey"]
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
