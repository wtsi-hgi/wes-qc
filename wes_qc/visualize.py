from typing import Optional

import hail as hl
import bokeh


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
