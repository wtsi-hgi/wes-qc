#get medians for number of variants per sample at each rf bin and for each cq
import hail as hl
import pyspark
import argparse
import matplotlib.pyplot as plt
import os
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


def get_vars_per_sample_per_bin_cq(mtfile: str, bins: list, consequences: list, plot_dir: str):
    '''
    get median numbers of variants per sample for any given consequence and random forest bin
    :param str mtfile: random forest and consequence annotated mtfile
    :param list bins: list of maximum bin numbers
    :param list consequences: vep consequences
    :param str plot_dir: output directory for plots
    '''
    mt = hl.read_matrix_table(mtfile)
    #filter to european samples
    mt_tmp = mt.filter_cols(mt.assigned_pop == 'EUR')
    sample_list = mt_tmp.cols().s.collect()
    samples = dict.fromkeys(sample_list)
    for s in samples:
        samples[s] = {}
    #create mts for indels and snps
    mt_snp = mt_tmp.filter_rows(hl.is_snp(mt_tmp.alleles[0], mt_tmp.alleles[1]))
    mt_indel = mt_tmp.filter_rows(hl.is_indel(mt_tmp.alleles[0], mt_tmp.alleles[1]))

    for consequence in consequences:
        print("Processing consequence " + consequence)
        mt_cq_snp_all = mt_snp.filter_rows(mt_snp.info.consequence == consequence)
        mt_cq_indel_all = mt_indel.filter_rows(mt_indel.info.consequence == consequence)

        for bin in bins:
            mt_cq_snp = mt_cq_snp_all.filter_rows(mt_cq_snp_all.info.rf_bin <= bin)
            mt_cq_indel = mt_cq_indel_all.filter_rows(mt_cq_indel_all.info.rf_bin <= bin)
            #run sample qc
            mt_cq_snp = hl.sample_qc(mt_cq_snp)
            mt_cq_indel = hl.sample_qc(mt_cq_indel)
            #extract into table
            sampleqc_ht_snp = mt_cq_snp.cols()
            sampleqc_ht_indel = mt_cq_indel.cols()
            #extract n_non_ref and sample ids - sample ids are extracted as we can't be sure that the sample qc 
            #table order will be the same each time
            non_ref_snp = sampleqc_ht_snp.sample_qc.n_non_ref.collect()
            non_ref_indel = sampleqc_ht_indel.sample_qc.n_non_ref.collect()
            samples_snp = sampleqc_ht_snp.s.collect()
            samples_indel = sampleqc_ht_indel.s.collect()
            snp_sample_counts = {samples_snp[i]: non_ref_snp[i] for i in range(len(samples_snp))}
            indel_sample_counts = {samples_indel[i]: non_ref_indel[i] for i in range(len(samples_indel))}

            for s_s in snp_sample_counts.keys():
                if not consequence in samples[s_s].keys():
                    samples[s_s][consequence] = {}
                if not bin in samples[s_s][consequence].keys():
                    samples[s_s][consequence][bin] = {}
                samples[s_s][consequence][bin]['snp'] = snp_sample_counts[s_s]

            for s_i in indel_sample_counts.keys():
                if not consequence in samples[s_i].keys():
                    samples[s_i][consequence] = {}
                if not bin in samples[s_i][consequence].keys():
                    samples[s_i][consequence][bin] = {}
                samples[s_i][consequence][bin]['indel'] = indel_sample_counts[s_s]

    #print to output file
    outfile = plot_dir + "/counts_per_sample.txt"
    header = ['sample']
    for cq in consequences:
        for bin in bins:
            header.append(cq + '_bin_' + str(bin) + '_snp')
            header.append(cq + '_bin_' + str(bin) + '_indel')
    with open(outfile, 'w') as o:
        o.write(('\t').join(header))
        o.write('\n')
        for s in samples.keys():
            outdata = []
            outdata.append(s)
            for cq in consequences:
                for bin in bins:
                    snp_count = samples[s][cq][bin]['snp']
                    indel_count = samples[s][cq][bin]['indel']
                    outdata.append(str(snp_count))
                    outdata.append(str(indel_count))
            o.write(("\t").join(outdata))
            o.write("\n")


def main():
    # set up
    args = get_options()
    inputs = parse_config()
    mtdir = inputs['matrixtables_lustre_dir']
    root_plot_dir = inputs['plots_dir_local']

    # initialise hail
    #tmp_dir = "hdfs://spark-master:9820/"
    tmp_dir = "file:///lustre/scratch123/qc/tmp"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    mtfile = mtdir + "mt_varqc_splitmulti_with_cq_and_rf_scores_" + args.runhash + ".mt"
    bins = list(range(1,101))
    consequences = ['missense_variant', 'synonymous_variant', 'frameshift_variant', 'inframe_deletion', 
    'inframe_insertion', 'splice_acceptor_variant', 'splice_donor_variant', 'stop_gained']

    plot_dir = root_plot_dir + "/variants_per_cq/" + args.runhash
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    get_vars_per_sample_per_bin_cq(mtfile, bins, consequences, plot_dir)


if __name__ == '__main__':
    main()
