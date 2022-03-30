# generate trusth sets for variant QC random forest
import hail as hl
import pyspark
import argparse
from wes_qc.utils.utils import parse_config

def get_options():
    '''
    Get options from the command line
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("--truth",
        help="create truth set HT", action="store_true")
    args = parser.parse_args()

    return args


def get_truth_ht(omni, mills, thousand_genomes, hapmap, truth_ht_file):
    '''
    Generate truth set HT
    :param str omni: file containing kgp_omni 1000 Genomes intersection Onni 2.5M array
    :param str mills: file containing Mills & Devine indels
    :param str thousand_genomes: file containing high confidence sites in 1000 genonmes
    :param str hapmap: file containing hapmap
    :param str truth_ht_file: output file name
    '''
    omni_ht = hl.read_table(omni)
    mills_ht = hl.read_table(mills)
    thousand_genomes_ht = hl.read_table(thousand_genomes)
    hapmap_ht = hl.read_table(hapmap)

    truth_ht = hapmap_ht.select(hapmap=True).join(
        omni_ht.select(omni=True), how="outer").join(
            thousand_genomes_ht.select(kgp_phase1_hc=True), how="outer").join(
                mills_ht.select(mills=True), how="outer").repartition(200, shuffle=False)
    truth_ht.write(truth_ht_file, overwrite=True)
    

def main():
    #set up
    args = get_options()
    inputs = parse_config()
    mtdir = inputs['matrixtables_lustre_dir']
    resourcedir = inputs['resource_dir']
    training_sets_dir = resourcedir + "training_sets/"

    # get truth set ht
    truth_ht_file = resourcedir + "truthset_table.ht"
    if args.truth:
        omni = training_sets_dir + "1000G_omni2.5.hg38.ht"
        mills = training_sets_dir + "Mills_and_1000G_gold_standard.indels.hg38.ht"
        thousand_genomes = training_sets_dir + "1000G_phase1.snps.high_confidence.hg38.ht"
        hapmap = training_sets_dir + "hapmap_3.3.hg38.ht"
        get_truth_ht(omni, mills, thousand_genomes, hapmap, truth_ht_file)

if __name__ == '__main__':
    main()
