#plots of variants per consequence or group of consequences
import argparse
import matplotlib.pyplot as plt
import pandas as pd
from wes_qc.utils.utils import parse_config

def get_options():
    '''
    Get options from the command line
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("--runhash", help="RF training run hash")
    parser.add_argument("--snv_bin", help="SNV cut off bin")
    parser.add_argument("--indel_bin", help="Indel cut off bin")
    args = parser.parse_args()
    if not args.runhash and args.snv_bin and args.indel_bin:
        print("--runhash, --snv_bin and --indel_bin_bin must be specified")
        exit(1)

    return args


def create_line_plots(plot_dir: str, medians_all: pd.Series, bins: list, consequences: dict):
    '''
    #create line plots for cumulative numbers or variants per sample for various consequences
    :param str plot_dir: output directory for plots
    :param pd.Series medians_all: medians for each cq/bin
    :param list bins: bins for which counts are available
    :param dict consequences: consequences/groups to plot with snp/indel
    '''
    for cq in consequences.keys():
        var_type = consequences[cq]
        medians_to_plot = []
        for b in bins:
            b = str(b)
            cq_col = cq + "_bin_" + b + "_" + var_type
            medians_to_plot.append(medians_all[cq_col])

        plt.plot(bins, medians_to_plot)
        plt.title("Median number of " +  cq + " variants per sample in each RF bin")
        plt.xlabel('Bin')
        plt.ylabel('Median variants per sample')
        plotfile = plot_dir + cq + "_cumulative_per_bin.png"
        plt.savefig(plotfile)


def create_violin_plots(plot_dir: str, df: pd.DataFrame, consequences: dict, snv_bin: str, indel_bin: str):
    '''
    create violin plots at specified cut off values
    :param str plot_dir: output directory for plots
    :param pd.DataFrame df: dataframe with counts per sample
    :param dict consequences: consequences/groups to plot with snp/indel
    :param str snv_bin: SNV bin cut off
    :param str indel_bin: indel bin cut off
    '''
    for cq in consequences:
        var_type = consequences[cq]
        if var_type == 'snp':
            plot_bin = snv_bin
        elif var_type == 'indel':
            plot_bin = indel_bin

        field_to_plot = cq + "_bin_" + plot_bin + "_" + var_type
        fig, ax = plt.subplots()
        ax.set_title("Variants per sample: " + cq + " at bin " + plot_bin)
        ax.set_ylabel("Variants per sample")
        violin_parts = ax.violinplot([ df[field_to_plot]], showmeans=False, showmedians=True, showextrema=False)
        for pc in violin_parts['bodies']:
            pc.set_facecolor('red')
            pc.set_edgecolor('black')

        plotfile = plot_dir + cq + "_violin.png"
        plt.savefig(plotfile)


def create_plots(plot_dir: str, bins: list, consequences: dict, snv_bin: str, indel_bin: str):
    '''
    create line plots of cumulative numbers of variants per bin per sample for various consequences
    :param str plot_dir: output directory for plots
    :param list bins: bins for which counts are available
    :param dict consequences: consequences/groups to plot with snp/indel
    :param str snv_bin: SNV bin cut off
    :param str indel_bin: indel bin cut off
    '''
    datafile = plot_dir + "/counts_per_sample.txt_test"
    df = pd.read_table(datafile)
    print(df.head())
    exit(0)
    #add additional fields to df  
    for b in bins:
        b = str(b)
        inframe_indel_bin = 'inframe_indels_bin_' + b + '_indel'
        inframe_ins_bin = 'inframe_insertion_bin_' + b + '_indel'
        inframe_del_bin = 'inframe_deletion_bin_' + b + '_indel'
        df[inframe_indel_bin] = df[inframe_ins_bin] + df[inframe_del_bin]
        lof_snv_bin = 'lof_snv_bin_' + b + '_snp'
        stop_gain_bin = 'stop_gained_bin_' + b + '_snp'
        splic_don_bin = 'splice_donor_variant_bin_' + b + '_snp'
        splic_acc_bin = 'splice_acceptor_variant_bin_' + b + '_snp'
        df[lof_snv_bin] = df[stop_gain_bin] + df[splic_acc_bin] + df[splic_don_bin]
    #get medians
    medians_all = df.median()

    #line plots
    create_line_plots(plot_dir, medians_all, bins, consequences)
    #violin_plots
    create_violin_plots(plot_dir, df, consequences, snv_bin, indel_bin)


def main():
    # set up
    args = get_options()
    inputs = parse_config()
    root_plot_dir = inputs['plots_dir_local']

    plot_dir = root_plot_dir + "/variants_per_cq/" + args.runhash + "/"

    bins = list(range(1,101))
    consequences = {
        'synonymous_variant': 'snp', 
        'missense_variant': 'snp', 
        'inframe_insertion': 'indel', 
        'inframe_deletion': 'indel', 
        'inframe_indels': 'indel', 
        'lof_snv': 'snp', 
        'frameshift_variant': 'indel'}
    create_plots(plot_dir, bins, consequences, args.snv_bin, args.indel_bin)


if __name__ == '__main__':
    main()