# create plots from binned random forest output
import hail as hl
import pyspark
import os
import argparse
import pandas as pd
from typing import Union, Dict, List, Set, Tuple
from bokeh.models import Plot, Row, Span
from bokeh.plotting import output_file, save
from gnomad.utils.plotting import *
from hail.plot import show, output_notebook
from wes_qc.utils.utils import parse_config


def get_options():
    '''
    Get options from the command line
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("--runhash", help="RF training run hash")
    args = parser.parse_args()
    if not args.runhash:
        print("--runhash must be specified")
        exit(1)

    return args


def set_plots_defaults(p: Plot, qc_plots_settings: dict) -> None:
    p.legend.label_text_font_size = qc_plots_settings['label_text_font_size']
    p.title.text_font_size = qc_plots_settings['title.text_font_size']
    p.axis.axis_label_text_font_size = qc_plots_settings['axis.axis_label_text_font_size']
    p.axis.axis_label_text_font_style = qc_plots_settings['axis.axis_label_text_font_style']
    p.axis.major_label_text_font_size = qc_plots_settings['axis.major_label_text_font_size']

    
def get_point_size_col(data: pd.Series, size_prop: str, qc_plots_settings: dict) -> pd.Series:
    """
    Given a data Series, returns the corresponding point size either:
    - Constant to qc_plots_settings['mean_point_size'] if `size_prop` is None
    - Radius proportional to data, if `size_prop` is 'radius'
    - Area proportional to data, if `size_prop` is 'area'

    Mean, min and max point  sizes are extracted from qc_plots_settings

    :param Series data: Input data series
    :param str size_prop: One of None, 'radius' or 'area'
    :return: Series with corresponding point size for each data point
    :rtype: Series
    """
    if size_prop is None:
        return pd.Series(len(data) * [qc_plots_settings['mean_point_size']])
    else:
        mean_data = np.mean(data)
        if size_prop == 'radius':
            return data.apply(lambda x: max(qc_plots_settings['min_point_size'], min(qc_plots_settings['max_point_size'], qc_plots_settings['mean_point_size'] * (x / mean_data))))
        elif size_prop == 'area':
            return data.apply(lambda x: max(qc_plots_settings['min_point_size'], min(qc_plots_settings['max_point_size'], qc_plots_settings['mean_point_size'] * np.pi * (np.sqrt(x / mean_data) / np.pi))))
        else:
            raise ValueError(f"{size_prop} is not a supported value for argument `size_prop`")


def plot_metric(df: pd.DataFrame,
                y_name: str,
                cols: List[str],
                qc_plots_settings: dict,
                y_fun: Callable[[pd.Series], Union[float, int]] = lambda x: x,
                cut: int = None,
                plot_all: bool = True,
                plot_bi_allelics: bool = True,
                plot_singletons: bool = True,
                plot_bi_allelic_singletons: bool = True,
                plot_adj: bool = False,
                colors: Dict[str, str] = None,
                link_cumul_y: bool = True,
                size_prop: str = 'area'
                ) -> Tabs:
    """
    Generic function for generating QC metric plots using a plotting-ready DataFrame (obtained from `get_binned_models_pd`)
    DataFrame needs to have a `rank_id` column, a `bin` column and a `model` column (contains the model name and needs to be added to binned table(s))

    This function generates scatter plots with the metric bin on x-axis and a user-defined function on the y-axis.
    The data for the y-axis function needs to from the columns specified in `cols`. The function is specified with the `y_fun` argument and data columns are access as a list.
    As an example, plotting Transition to transversion ratio is done as follows::

        plot_metric(snvs, 'Ti/Tv', ['n_ti', 'n_tv'], y_fun=lambda x: x[0]/x[1], colors=colors)

    In this command, `x[0]` correspond to the  first column selected (`'n_ti'`)  and `x[1]` to the second (`'n_tv'`).


    This function plots a tab for each of the plot condition(s) selected: all, bi-allelics, bi-allelic singletons.
    Within each tab, each row contains a non-cumulative and a cumulative plot of the bins / values.
    If `plot_adj` is set, then an extra row is added plotting only variants in release samples where AC_ADJ>0. The bin for these sites is computed based on those variants only.

    :param pd.DataFrame df: Input data
    :param str y_name: Name of the metric plotted on the y-axis
    :param list of str cols: Columns used to compute the metric plotted
    :pararm dict qc_plots_settings: Plot settings
    :param callable y_fun: Function to apply to the columns to generate the metric
    :param int cut: Where to draw the bin cut
    :param bool plot_all: Whether to plot a tab with all variants
    :param bool plot_bi_allelics: Whether to plot a tab with bi-allelic variants only
    :param bool plot_singletons: Whether to plot a tab with singleton variants only
    :param bool plot_bi_allelic_singletons:  Whether to plot a tab with bi-allelic singleton variants only
    :param bool plot_adj: Whether to plot additional rows with adj variants in release samples only
    :param dict of str -> str colors: Mapping of model name -> color
    :param bool link_cumul_y: If set, y-axes of cumulative and non-cumulative plots are linked
    :param str size_prop: Either 'size' or 'area' can be specified. If either is specified, the points will be sized proportionally to the amount of data in that point.
    :return: Plot
    :rtype: Tabs
    """
    def get_row(df: pd.DataFrame, y_name: str, cols: List[str], y_fun: Callable[[pd.Series], Union[float, int]], titles: List[str], link_cumul_y: bool, cut: int = None) -> Row:
        """
        Generates a single row with two plots: a regular scatter plot and a cumulative one.
        Both plots have bins on the x-axis. The y-axis is computed by applying the function `y_fun` on the columns `cols`.

        Data source is shared between the two plots so that highlighting / selection is linked.
        X-axis is shared between the two plots.
        Y-axus is shared if `link_cumul_y` is `True`

        """

        def get_plot(data_source: ColumnDataSource, y_name: str, y_col_name: str, titles: List[str], data_ranges: Tuple[DataRange1d, DataRange1d], cut: int = None) -> Plot:
            """
            Generates a single scatter plot panel
            """

            p = figure(
                title=titles[0],
                x_axis_label='bin',
                y_axis_label=y_name,
                tools="save,pan,box_zoom,reset,wheel_zoom,box_select,lasso_select,help,hover")
            p.x_range = data_ranges[0]
            p.y_range = data_ranges[1]

            if cut:
                p.add_layout(Span(location=cut, dimension='height', line_color='red', line_dash='dashed'))

            # Add circles layouts one model at a time, so that no default legend is generated.
            # Because data is in the same ColumnDataSource, use a BooleanFilter to plot each model separately
            circles = []
            for model in set(data_source.data['model']):
                view = CDSView(source=data_source, filters=[BooleanFilter([x == model for x in data_source.data['model']])])
                circles.append((model, [p.circle('bin', y_col_name, color='_color', size='_size', source=data_source, view=view)]))

            p.select_one(HoverTool).tooltips = [('model', '@model'),
                                                ('bin', '@bin'),
                                                (y_name, f'@{y_col_name}'),
                                                ('min_score', '@min_score'),
                                                ('max_score', '@max_score'),
                                                ('n_data_points', '@_n')
                                                ] + [(col, f'@{col}') for col in cols]
            set_plots_defaults(p, qc_plots_settings)

            # Add legend above the plot area
            legend = Legend(items=circles, orientation='horizontal', location=(0, 0), click_policy="hide")
            p.add_layout(legend, 'above')

            # Add subtitles if any
            for title in titles[1:]:
                p.add_layout(Title(text=title, text_font_size=qc_plots_settings['subtitle.text_font_size']), 'above')

            return p
        # Compute non-cumulative values by applying `y_fun`
        df['non_cumul'] = df[cols].apply(y_fun, axis=1)

        # Compute cumulative values for each of the data columns
        for col in cols:
            df[f'{col}_cumul'] = df.groupby('model').aggregate(np.cumsum)[col]
        df['cumul'] = df[[f'{col}_cumul' for col in cols]].apply(y_fun, axis=1)

        # Create data ranges that are either shared or distinct depending on the y_cumul parameter
        non_cumul_data_ranges = (DataRange1d(), DataRange1d())
        cumul_data_ranges = non_cumul_data_ranges if link_cumul_y else (non_cumul_data_ranges[0], DataRange1d())
        data_source = ColumnDataSource(df)

        return Row(get_plot(data_source, y_name, 'non_cumul', titles, non_cumul_data_ranges, cut),
                   get_plot(data_source, y_name, 'cumul', [titles[0] + ', cumulative'] + titles[1:], cumul_data_ranges, cut))

    
    def prepare_pd(df: pd.DataFrame, cols: List[str], colors: Dict[str, str] = {}, size_prop: str = None):
        """
        Groups a pandas DataFrame by model and bin while keeping relevant columns only.
        Adds 3 columns used for plotting:
        1. A _color column column
        2. A _n column containing the number of data points
        3. A _size column containing the size of data points based on the `size_prop` and `qc_plot_settings` parameters
        """
        df = df.groupby(['model', 'bin']).agg({**{col: np.sum for col in cols},
                                               'min_score': np.min, 'max_score': np.max})
        df = df.reset_index()
        df['_color'] = [colors.get(x, 'gray') for x in df['model']]
        df['_n'] = np.sum(df[cols], axis=1)
        df['_size'] = get_point_size_col(df['_n'], size_prop, qc_plots_settings)
        return df

    colors = colors if colors is not None else {}
    tabs = []
    adj_strats = ['', 'adj_'] if plot_adj else ['']

    if plot_all:
        children = []
        for adj in adj_strats:
            titles = [y_name, 'Adj variants (adj rank)' if adj else 'All variants']
            plot_df = prepare_pd(df.loc[df.rank_id == f'{adj}rank'], cols, colors, size_prop)
            if len(plot_df) > 0:
                children.append(get_row(plot_df, y_name, cols, y_fun, titles, link_cumul_y, cut))
            else:
                print('No data found for plot: {}'.format('\t'.join(titles)))

        if children:
            tabs.append(Panel(child=Column(children=children), title='All'))

   
    if plot_singletons:
        children = []
        for adj in adj_strats:
            for singleton_rank in ['', 'singleton_']:
                titles = [y_name, 'Singletons ({} rank)'.format('overall' if not adj and not singleton_rank else " ".join([adj[:-1], singleton_rank[:-1]]).lstrip())]
                plot_df = prepare_pd(df.loc[df.singleton & (df.rank_id == f'{adj}{singleton_rank}rank')], cols, colors, size_prop)
                if len(plot_df) > 0:
                    children.append(get_row(plot_df, y_name, cols, y_fun, titles, link_cumul_y, cut))
                else:
                    print('No data found for plot: {}'.format('\t'.join(titles)))

        if children:
            tabs.append(Panel(child=Column(children=children), title='Singletons'))

    return Tabs(tabs=tabs)


def create_plots(bin_htfile: str, plot_dir: str, run_hash: str, qc_plots_settings: dict):
    '''
    Create variant QC plots
    :param str bin_htfile: Hail table file containing binned ranfom forest output
    :param str plot_dir: Output plot directory
    :param str run_hash: hash of random forest run
    :param dict qc_plots_settings: qc plots settings
    '''
    #convert ht to pandas df
    ht = hl.read_table(bin_htfile)
    #drop de_novo_high_quality as this field is sometimes NA and messes with conversion to pands
    ht = ht.key_by('rank_id', 'contig', 'snv', 'bi_allelic', 'singleton', 'trans_singletons', 'de_novo_medium_quality', 'de_novo_synonymous', 'bin')
    ht = ht.drop('de_novo_high_quality')

    ht = ht.group_by(
            *[k for k in ht.key if k != 'contig']
        ).aggregate(
            min_score=hl.agg.min(ht.min_score),
            max_score=hl.agg.max(ht.max_score),
            **{x: hl.agg.sum(ht[x]) for x in ht.row_value if x not in ['min_score', 'max_score']}
        )
    ht = ht.annotate(model=run_hash)
    df = ht.to_pandas()
    #indel and SNVs dataframes
    snvs = df[df['snv']==True]
    indels = df[df['snv']==False]
    #plots
    model = run_hash
    colors = {model:'blue'}
    #plot transmitted singletons
    plotfile = plot_dir + "transmitted_singletons.html"
    tabs = plot_metric(snvs, 'n_trans_singletons', ['n_trans_singletons_synonymous_algorithm'], qc_plots_settings, y_fun=lambda x: x[0], plot_bi_allelics=False, plot_singletons=False, plot_bi_allelic_singletons=False, colors=colors)
    output_file(filename=plotfile)
    save(tabs)
    #plot untransmitted singletons
    plotfile = plot_dir + "untransmitted_singletons.html"
    tabs = plot_metric(snvs, 'n_untrans_singletons', ['n_untrans_singletons_synonymous_algorithm'], qc_plots_settings, y_fun=lambda x: x[0], plot_bi_allelics=False, plot_singletons=False, plot_bi_allelic_singletons=False, colors=colors)
    output_file(filename=plotfile)
    save(tabs)
    #plot transmitted/untransmitted ratio
    plotfile = plot_dir + "transmitted_untransmitted.html"
    tabs = plot_metric(snvs, 'trans_untrans_ratio', ['n_trans_singletons_synonymous_algorithm', 'n_untrans_singletons_synonymous_algorithm'], qc_plots_settings, y_fun=lambda x: x[0]/x[1], plot_bi_allelics=False, plot_singletons=False, plot_bi_allelic_singletons=False, colors=colors)
    output_file(filename=plotfile)
    save(tabs)
    #plot number of insertions
    plotfile = plot_dir + "n_insertions.html"
    tabs = plot_metric(indels, 'n_ins', ['n_ins'], qc_plots_settings, y_fun=lambda x: x[0], plot_bi_allelics=False, plot_singletons=False, plot_bi_allelic_singletons=False, colors=colors)
    output_file(filename=plotfile)
    save(tabs)
    #plot number of deletions
    plotfile = plot_dir + "n_deletions.html"
    tabs = plot_metric(indels, 'n_del', ['n_del'], qc_plots_settings, y_fun=lambda x: x[0], plot_bi_allelics=False, plot_singletons=False, plot_bi_allelic_singletons=False, colors=colors)
    output_file(filename=plotfile)
    save(tabs)
    #plot Ti/Tv ratio
    plotfile = plot_dir + "r_ti_tv.html"
    tabs = plot_metric(snvs, 'Ti/Tv ratio', ['n_ti', 'n_tv'], qc_plots_settings, y_fun=lambda x: x[0]/x[1], plot_bi_allelics=False, plot_singletons=False, plot_bi_allelic_singletons=False, colors=colors)
    output_file(filename=plotfile)
    save(tabs)
    #plot 1kg high confidence SNVs
    plotfile = plot_dir + "kg_snv.html"
    tabs = plot_metric(snvs, 'n_kgp_high_conf_snvs', ['n_kgp_high_conf_snvs'], qc_plots_settings, y_fun=lambda x: x[0], plot_bi_allelics=False, plot_singletons=False, plot_bi_allelic_singletons=False, colors=colors)
    output_file(filename=plotfile)
    save(tabs)
    #plot omni SNVs
    plotfile = plot_dir + "omni_snv.html"
    tabs = plot_metric(snvs, 'n_omni', ['n_omni'], qc_plots_settings, y_fun=lambda x: x[0], plot_bi_allelics=False, plot_singletons=False, plot_bi_allelic_singletons=False, colors=colors)
    output_file(filename=plotfile)
    save(tabs)
    #plot mills indels
    plotfile = plot_dir + "mills_indels.html"
    tabs = plot_metric(indels, 'n_mills', ['n_mills'], qc_plots_settings, y_fun=lambda x: x[0], plot_bi_allelics=False, plot_singletons=False, plot_bi_allelic_singletons=False, colors=colors)
    output_file(filename=plotfile)
    save(tabs)
    #plot hapmap snvs
    plotfile = plot_dir + "hapmap_snvs.html"
    tabs = plot_metric(snvs, 'n_hapmap_snvs', ['n_hapmap'], qc_plots_settings, y_fun=lambda x: x[0], plot_bi_allelics=False, plot_singletons=False, plot_bi_allelic_singletons=False, colors=colors)
    output_file(filename=plotfile)
    save(tabs)
    #plot hapmap indels
    plotfile = plot_dir + "hapmap_indels.html"
    tabs = plot_metric(indels, 'n_hapmap_indels', ['n_hapmap'], qc_plots_settings, y_fun=lambda x: x[0], plot_bi_allelics=False, plot_singletons=False, plot_bi_allelic_singletons=False, colors=colors)
    output_file(filename=plotfile)
    save(tabs)
    #plot fail hard filters snvs
    plotfile = plot_dir + "fail_hard_filters_snvs.html"
    tabs = plot_metric(snvs, 'fail_hard_filters_snvs', ['fail_hard_filters'], qc_plots_settings, y_fun=lambda x: x[0], plot_bi_allelics=False, plot_singletons=False, plot_bi_allelic_singletons=False, colors=colors)
    output_file(filename=plotfile)
    save(tabs)
    #plot fail hard filters indels
    plotfile = plot_dir + "fail_hard_filters_indels.html"
    tabs = plot_metric(snvs, 'fail_hard_filters_indels', ['fail_hard_filters'], qc_plots_settings, y_fun=lambda x: x[0], plot_bi_allelics=False, plot_singletons=False, plot_bi_allelic_singletons=False, colors=colors)
    output_file(filename=plotfile)
    save(tabs)


def main():
    # set up
    args = get_options()
    inputs = parse_config()
    rf_dir = inputs['var_qc_rf_dir']
    root_plot_dir = inputs['plots_dir_local']

    qc_plots_settings = {
    'mean_point_size': 4.0,
    'min_point_size': 1.0,
    'max_point_size': 16.0,
    'label_text_font_size': "14pt",
    'title.text_font_size': "16pt",
    'subtitle.text_font_size': "14pt",
    'axis.axis_label_text_font_size': "16pt",
    'axis.axis_label_text_font_style': "normal",
    'axis.major_label_text_font_size': "14pt"
    }
    plot_dir = root_plot_dir + "variant_qc/" + args.runhash + "/"
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    bin_htfile = rf_dir + args.runhash + "/_rf_result_ranked_BINS.ht"
    create_plots(bin_htfile, plot_dir, args.runhash, qc_plots_settings)


if __name__ == '__main__':
    main()