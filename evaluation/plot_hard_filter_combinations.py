#create interactive plots showing hard filter combinations
import plotly.express as px
import pandas as pd
from wes_qc.utils.utils import parse_config


def make_plot(df: pd.DataFrame, x: str, y:str, outfile: str):
    '''
    Creata plot from apanda df, x column, y column
    :param pd.dDtaFrame df: panda dataframe
    :param str x: column for x axis
    :param str y: column for y axis
    :param str outfile: plot output file
    '''
    fig = px.scatter(df, x=x, y=y, color = "filter")
    fig.write_html(outfile)

def create_plots(snv_results_file: str, indel_results_file: str, outdir: str):
    '''
    Create plots
    :param str snv_results_file: snv results file path
    :param str indel_results_file: indel results file path
    :param str outdir: output directory for plots
    '''
    snv_df = pd.read_csv(snv_results_file, sep = "\t")
    indel_df = pd.read_csv(indel_results_file, sep = "\t")
    
    snv_tp_fp_file = outdir + "snv_tp_fp.html"
    make_plot(snv_df, "TP", "FP", snv_tp_fp_file)


def main():
    # set up
    inputs = parse_config()
    plot_dir = inputs['plots_dir_local']

    snv_results_file = plot_dir + "genotype_hard_filter_comparison_snv_v1.txt"
    indel_results_file = plot_dir + "genotype_hard_filter_comparison_indel_v1.txt"
    outdir = plot_dir + "hard_filter_evaluation/"

    create_plots(snv_results_file, indel_results_file, outdir)

if __name__ == '__main__':
    main()