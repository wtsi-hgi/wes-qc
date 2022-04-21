#take file with the best hit for each sample from gtcheck
#for those with a good score, does their best hit match the expected id from the metadata file?
#are there any ids in the genotyping data that are the best hit for more than one WES sample?
#identify those with no good hit (score < 0.05). What genotyping sample do we expect them to correspond to 
#and is that sample in the gneotyping file?


def find_plink_samples(plinksamplesfile):
    '''
    create a list of samples that are in the plink data and parse to remove repetition (samples in plink are in format 123A_123A)
    '''
    samples = []
    with open(plinksamplesfile, 'r') as f:
        samples = f.readlines()
        samples = [x.strip() for x in samples]
        samples = [x.split("_")[0] for x in samples]

    return samples


def main():
    gtcheck_output_file = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/whole_exome_output/gtcheck_best_hits.txt'
    plinksamplesfile = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/plink_samples.txt'
    metadata_file = '/lustre/scratch119/humgen/projects/birth_cohort_wes/genestack/all_samples_with_proceed_and_seq_info_and_warehouse_info.txt'
    plink_samples = find_plink_samples(plinksamplesfile)
    print(plink_samples)


if __name__ == '__main__':
    main()