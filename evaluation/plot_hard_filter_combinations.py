#create interactive plots showing hard filter combinations
import plotly.express as px
import pandas as pd
from wes_qc.utils.utils import parse_config


def make_plot(df: pd.DataFrame, x: str, y:str, outfile: str, vartype: str):
    '''
    Creata plot from apanda df, x column, y column
    :param pd.dDtaFrame df: panda dataframe
    :param str x: column for x axis
    :param str y: column for y axis
    :param str outfile: plot output file
    :param str vartype: Variant type
    '''
    symbolmap = {'5': 'circle', '10': 'cross'}
    if x == "t_u_ratio":
        plottitle = (" ").join([vartype, "transmitted/untransmitted synonymous singletons", y])
    else:
        plottitle = (" ").join([vartype, x, y])
    fig = px.scatter(df, x=x, y=y, color = "bin", symbol = 'DP', symbol_map = symbolmap, hover_data = ['GQ', 'AB'], title = plottitle)
    fig.write_html(outfile)

def create_plots(snv_results_file: str, indel_results_file: str, outdir: str):
    '''
    Create plots
    :param str snv_results_file: snv results file path
    :param str indel_results_file: indel results file path
    :param str outdir: output directory for plots
    '''
    snv_df = pd.read_csv(snv_results_file, sep = "\t")
    snv_df[['w', 'bin', 'x', 'DP', 'y', 'GQ', 'z', 'AB']] = snv_df['filter'].str.split('_', expand=True)
    snv_df.drop(['w', 'x', 'y', 'z'], axis = 1, inplace = True)

    indel_df = pd.read_csv(indel_results_file, sep = "\t")
    indel_df[['w', 'bin', 'x', 'DP', 'y', 'GQ', 'z', 'AB']] = indel_df['filter'].str.split('_', expand=True)
    indel_df.drop(['w', 'x', 'y', 'z'], axis = 1, inplace = True)
    
    snv_tp_fp_file = outdir + "snv_tp_fp.html"
    make_plot(snv_df, "TP", "FP", snv_tp_fp_file, 'SNV')

    snv_tu_tp_file = outdir + "snv_tu_tp.html"
    make_plot(snv_df, "t_u_ratio", "TP", snv_tu_tp_file, 'SNV')

    snv_prec_recall_file = outdir + "snv_prec_recall.html"
    make_plot(snv_df, "recall", "precision", snv_prec_recall_file, 'SNV')

    indel_tp_fp_file = outdir + "indel_tp_fp.html"
    make_plot(indel_df, "TP", "FP", indel_tp_fp_file, 'Indel')

    indel_prec_recall_file = outdir + "indel_prec_recall.html"
    make_plot(indel_df, "recall", "precision", indel_prec_recall_file, 'Indel')

    indel_prec_recall_file_frameshift = outdir + "indel_prec_recall_frameshift.html"
    make_plot(indel_df, "recall_frameshift", "precision_frameshift", indel_prec_recall_file_frameshift, 'Indel')

    indel_prec_recall_file_inframe = outdir + "indel_prec_recall_inframe.html"
    make_plot(indel_df, "recall_inframe", "precision_inframe", indel_prec_recall_file_inframe, 'Indel')


def main():
    # set up
    inputs = parse_config()
    plot_dir = inputs['plots_dir_local']

    # snv_results_file = plot_dir + "genotype_hard_filter_comparison_snv_v1.txt"
    # indel_results_file = plot_dir + "genotype_hard_filter_comparison_indel_v1.txt"
    # snv_results_file = plot_dir + "genotype_hard_filter_comparison_snv.txt"
    # indel_results_file = plot_dir + "genotype_hard_filter_comparison_indel.txt"
    snv_results_file = plot_dir + "86f80184_genotype_hard_filter_comparison_snv_fewer_combs_v2.txt"
    indel_results_file = plot_dir + "86f80184_genotype_hard_filter_comparison_indel_fewer_combs_v2.txt"
    outdir = plot_dir + "hard_filter_evaluation/230329/"

    create_plots(snv_results_file, indel_results_file, outdir)

if __name__ == '__main__':
    main()