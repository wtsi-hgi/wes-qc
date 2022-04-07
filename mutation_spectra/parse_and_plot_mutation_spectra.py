#parse the bcftools stats output and create mutation specta for a list of samples
#also output a file with fraction of each substitution for each sample

import argparse
import os
import gzip

def get_options():
    '''
    gets options from the command line
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--sample_list", required=True,
                            help="sample list")
    parser.add_argument("-d", "--bcftools_stats_directory", required=True,
                            help="path to bcftools stats directory")
    parser.add_argument("-o", "--output_directory", required=True,
                            help="path to bcftools output directory")
    args = parser.parse_args()

    if not os.path.isdir(args.output_directory):
        print("Directory not found: " + args.output_directory + " creating new directory")
        os.makedirs(args.output_directory)

    if os.path.isfile(args.sample_list) and os.path.isdir(args.bcftools_stats_directory):
        return args
    else:
        if not os.path.isfile(args.sample_list):
            print("File not found: " + args.sample_list)
        if not os.path.isdir(args.bcftools_stats_directory):
            print("Directory not found: " + args.bcftools_stats_directory)
        exit(1)


def parse_sample_list(sample_list):
    '''
    parse sample list file and return a list of samples
    '''
    with open(sample_list, 'r') as f:
        samples = f.readlines()
        samples = [s.rstrip() for s in samples]

    return samples


def parse_bcftools_stats(samples, bcftoos_stats_dir, outdir):
    '''
    for each per-person, per-chromsome output file, combine counts of 
    substititions, work out fraction, create a summary and write this to an output file
    '''
    chroms = list(range(1,23))
    chroms = [str(x) for x in chroms]
    chroms.extend(['X', 'Y'])
    substitutions = ['A>C', 'A>G', 'A>T', 'C>A', 'C>G',
                     'C>T', 'G>A', 'G>C', 'G>T', 'T>A', 'T>C', 'T>G']
    props_per_sample = {}

    for s in samples:
        counts_per_sample = {'A>C': 0, 'A>G': 0, 'A>T': 0, 'C>A': 0, 'C>G': 0, 'C>T': 0,
                             'G>A': 0, 'G>C': 0, 'G>T': 0, 'T>A': 0, 'T>C': 0, 'T>G': 0, 'total': 0}        
        for c in chroms:
            statsfile = bcftoos_stats_dir + "/" + s + "_" + c + ".stats.gz"
            with gzip.open(statsfile, 'rt') as zf:
                for line in zf:
                    if line.startswith('ST'):
                        linedata = line.split()
                        subs = linedata[2]
                        count = int(linedata[3])
                        counts_per_sample[subs] += count
                        counts_per_sample['total'] += count
        props_per_sample[s] = {}
        for st in substitutions:
            props_per_sample[s][st] = counts_per_sample[st]/counts_per_sample['total']

    summaryfile = outdir + "/proportions_per_person.txt"
    with open(summaryfile, 'w') as o:
        header = sample + "\t" + ("\t").join(substitutions)
        o.write(header)
        o.write("\n")
        for sample in props_per_sample.keys():
            outdata = [sample]
            for st in substitutions:
                num = "{:.3f}".format(props_per_sample[sample][st])
                outdata.append(num)
            o.write(("\t").join(outdata))
            o.write("\n")

    return props_per_sample


def create_plots(props_per_sample, outdir):
    '''
    create a plot for each sample
    '''
    pass


def main():
    args = get_options()
    samples = parse_sample_list(args.sample_list)
    props_per_sample = parse_bcftools_stats(samples, args.bcftools_stats_directory, args.output_directory)
    create_plots(props_per_sample, args.output_directory)


if __name__ == "__main__":
    main()
