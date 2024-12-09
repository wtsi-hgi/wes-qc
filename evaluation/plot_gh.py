#create interactive plots showing hard filter combinations
import plotly.express as px
import pandas as pd
#from utils.utils import parse_config
import json

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
    fig = px.scatter(df, x=x, y=y, color = "bin", symbol = 'DP', symbol_map = symbolmap, hover_data = ['GQ', 'AB'], title = plottitle, facet_col='missing')
    #fig = px.scatter(df, x=x, y=y, color = "bin", symbol = 'DP', symbol_map = symbolmap, size='AB', hover_data = 'GQ', title = plottitle, facet_col='missing')
    fig.write_html(outfile)

def create_plots(json_file1: str, json_file2: str, outdir: str):
    '''
    Create plots
    :param str snv_results_file: snv results file path
    :param str indel_results_file: indel results file path
    :param str outdir: output directory for plots
    '''
    with open(json_file1, 'r') as file:
        json_data1 = json.load(file)
    indel_df=pd.DataFrame.from_dict(json_data1['indel'], orient='index').reset_index()
    indel_df.rename(columns={'index': 'filter'}, inplace=True)
    with open(json_file2, 'r') as file:
        json_data2 = json.load(file)
    snv_df=pd.DataFrame.from_dict(json_data2['snv'], orient='index').reset_index()
    snv_df.rename(columns={'index': 'filter'}, inplace=True)
    #snv_df = pd.read_json(snv_results_file, lines=True)
    snv_df[['v', 'bin', 'w', 'DP', 'x', 'GQ', 'y', 'AB', 'z', 'missing']] = snv_df['filter'].str.split('_', expand=True)
    snv_df.drop(['v', 'w', 'x', 'y', 'z'], axis = 1, inplace = True)

    #indel_df = pd.read_json(indel_results_file, lines=True)
    indel_df[['v', 'bin', 'w', 'DP', 'x', 'GQ', 'y', 'AB', 'z', 'missing']] = indel_df['filter'].str.split('_', expand=True)
    indel_df.drop(['v', 'w', 'x', 'y', 'z'], axis = 1, inplace = True)
    
    indel_df['TP'] = (indel_df['TP'] / json_data1['indel_total_tp']) * 100
    indel_df['FP'] = (indel_df['FP'] / json_data1['indel_total_fp']) * 100
    snv_df['TP'] = (snv_df['TP'] / json_data2['snv_total_tp']) * 100
    snv_df['FP'] = (snv_df['FP'] / json_data2['snv_total_fp']) * 100

    snvf=outdir+"/snv_df.tsv"
    indelf=outdir+"/indel_df.tsv"
    snv_df.to_csv(snvf, index=False, sep="\t")
    indel_df.to_csv(indelf, index=False, sep="\t")

    #max_snv = snv_df["precision"].max()
    #min_snv = snv_df["precision"].min()
    #max_indel = indel_df[["precision", "precision_frameshift", "precision_inframe"]].max(axis=1).max()
    #min_indel = indel_df[["precision", "precision_frameshift", "precision_inframe"]].min(axis=1).min()

    #min_precision = min(min_snv, min_indel)-0.01
    #max_precision = max(max_snv, max_indel)+0.01

    #max_snv = snv_df["recall"].max()
    #min_snv = snv_df["recall"].min()
    #max_indel = indel_df[["recall", "recall_frameshift", "recall_inframe"]].max(axis=1).max()
    #min_indel = indel_df[["recall", "recall_frameshift", "recall_inframe"]].min(axis=1).min()

    #min_recal = min(min_snv, min_indel)-0.01
    #max_recal = max(max_snv, max_indel)+0.01


    snv_tp_fp_file = outdir + "/snv_tp_fp.html"
    make_plot(snv_df, "TP", "FP", snv_tp_fp_file, 'SNV')

    snv_tu_tp_file = outdir + "/snv_tu_tp.html"
    make_plot(snv_df, "t_u_ratio", "TP", snv_tu_tp_file, 'SNV')

    snv_prec_recall_file = outdir + "/snv_prec_recall.html"
    make_plot(snv_df, "recall", "prec", snv_prec_recall_file, 'SNV')

    snv_tp_fp_file = outdir + "/snv_mendel_tp_missing.html"
    make_plot(snv_df, "mendel_errors", "TP", snv_tp_fp_file, 'SNV')

    indel_tp_fp_file = outdir + "/indel_tp_fp.html"
    make_plot(indel_df, "TP", "FP", indel_tp_fp_file, 'Indel')

    indel_prec_recall_file = outdir + "/indel_prec_recall.html"
    make_plot(indel_df, "recall", "prec", indel_prec_recall_file, 'Indel')

    indel_prec_recall_file_frameshift = outdir + "/indel_prec_recall_frameshift.html"
    make_plot(indel_df, "recall_frameshift", "prec_frameshift", indel_prec_recall_file_frameshift, 'Indel')

    indel_prec_recall_file_inframe = outdir + "/indel_prec_recall_inframe.html"
    make_plot(indel_df, "recall_inframe", "prec_inframe", indel_prec_recall_file_inframe, 'Indel')

    indel_tp_fp_file = outdir + "/indel_mendel_tp_missing.html"
    make_plot(indel_df, "mendel_errors", "TP", indel_tp_fp_file, 'Indel')


def main():
    # set up
    #inputs = parse_config()
    plot_dir = '/lustre/scratch123/mdt2/projects/gnh_industry/Genes_and_Health_2024_05_55k/qc/plots/'

    # snv_results_file = plot_dir + "genotype_hard_filter_comparison_snv_v1.txt"
    # indel_results_file = plot_dir + "genotype_hard_filter_comparison_indel_v1.txt"
    #snv_results_file = plot_dir + "genotype_hard_filter_comparison_snv.txt"
    #indel_results_file = plot_dir + "genotype_hard_filter_comparison_indel.txt"
    #snv_results_file = plot_dir + "genotype_hard_filter_comparison_snv_fewer_combs_v3.test3.txt"
    #indel_results_file = plot_dir + "genotype_hard_filter_comparison_indel_fewer_combs_v3.test3.txt"
    json_file1='/lustre/scratch123/mdt2/projects/gnh_industry/Genes_and_Health_2024_05_55k/qc/matrixtables/0d975fbd/evaluation.indel.json'
    json_file2='/lustre/scratch123/mdt2/projects/gnh_industry/Genes_and_Health_2024_05_55k/qc/matrixtables/0d975fbd/evaluation.snv.json'
    outdir = plot_dir + "new_plots"

    create_plots(json_file1, json_file2, outdir)

if __name__ == '__main__':
    main()

