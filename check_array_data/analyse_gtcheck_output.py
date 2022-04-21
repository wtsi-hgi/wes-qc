#take file with the best hit for each sample from gtcheck
#for those with a good score, does their best hit match the expected id from the metadata file?
#are there any ids in the genotyping data that are the best hit for more than one WES sample?
#identify those with no good hit (score < 0.05). What genotyping sample do we expect them to correspond to 
#and is that sample in the gneotyping file?


def find_plink_samples(plinksamples_file):
    '''
    create a list of samples that are in the plink data and parse to remove repetition (samples in plink are in format 123A_123A)
    '''
    samples = []
    with open(plinksamples_file, 'r') as f:
        samples = f.readlines()
        samples = [x.strip() for x in samples]
        samples = [x.split("_")[0] for x in samples]

    return samples


def parse_metadata(metadata_file):
    '''
    create a dict of WES sample id and genotyping sample id
    '''
    sample_map = {}
    with open(metadata_file, 'r') as f:
        lines = f.readlines()
        for l in lines:
            if l.startswith('sangersampleid'):
                continue
            ldata = l.split()
            if ldata[0].startswith('Z'):
                wes_sample = ldata[0]
            elif ldata[32].startswith('EGA'):
                wes_sample = ldata[32]
            else:
                print("No WES sample id found in line: " + l)
                exit(1)
            plink_sample = ldata[14]
            sample_map[wes_sample] = plink_sample
            print(sample_map)
            exit(0)


    return sample_map


def main():
    gtcheck_output_file = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/whole_exome_output/gtcheck_best_hits.txt'
    plinksamples_file = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/plink_samples.txt'
    metadata_file = '/lustre/scratch119/humgen/projects/birth_cohort_wes/genestack/all_samples_with_proceed_and_seq_info_and_warehouse_info.txt'

    plink_samples = find_plink_samples(plinksamples_file)
    sample_map = parse_metadata(metadata_file)



if __name__ == '__main__':
    main()