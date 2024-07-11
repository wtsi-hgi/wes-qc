import hail as hl
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns


def visualize_pca(pca_scores: hl.Table, plot_file: str, n_pca: int = 3) -> None:
    print("Plotting PCA components form Hail table")
    sns_data = pd.DataFrame.from_records(pca_scores.to_pandas().scores)[range(n_pca)]
    sns_data["s"] = pca_scores.s.collect()
    p = sns.pairplot(sns_data, vars=sns_data.columns[range(n_pca)], plot_kws={"s": 2})
    plt.savefig(plot_file)
