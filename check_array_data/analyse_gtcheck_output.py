# take file with the best hit for each sample from gtcheck
# for those with a good score, does their best hit match the expected id from the metadata file?
# are there any ids in the genotyping data that are the best hit for more than one WES sample?
# identify those with no good hit (score < 0.05). What genotyping sample do we expect them to correspond to
# and is that sample in the gneotyping file?


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


def parse_metadata(metadata_file, gtcheck_duplicates_file):
    '''
    create a dict of WES sample id and genotyping sample id
    identify genotyping sample ids that correspond to more than one WES sample and print to file
    '''
    sample_map = {}
    plink_to_ega = {}
    duplicates = {}
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
                
            plink_sample = ldata[14]
            sample_map[wes_sample] = plink_sample
            if not plink_sample in plink_to_ega.keys():
                plink_to_ega[plink_sample] = wes_sample
            else:
                if not plink_sample in duplicates.keys():
                    duplicates[plink_sample] = [plink_to_ega[plink_sample]]
                duplicates[plink_sample].append[wes_sample]

    with open(gtcheck_duplicates_file, 'w') as o:
        for s in duplicates.keys():
            o.write(s + (", ").join(duplicates[s]))
            o.write("\n")

    return sample_map


def parse_gtcheck_output(gtcheck_output_file, plink_samples, sample_map, gtcheck_dodgy_samples_file, gtcheck_mismatches_file):
    '''
    parse gtcheck output
    for those with a good hit, check that the hit was with the expected id
    for those with no good hit check if we expect the sample to be present in the plink file
    check for families in plink that appear to be duplicated in WES
    '''
    mismatches = []
    dodgy_samples = {}
    with open(gtcheck_output_file, 'r') as g:
        lines = g.readlines()
        for l in lines:
            if l.startwith('#'):
                continue
            ldata = l.split()
            wes_sample = ldata[0]
            plink_sample = ldata[1].split("_")[0]
            score = float(ldata[4])
            # score of < 0.05 seen as good so check if the samples match - if they don't it is a mismatch
            if score < 0.05:
                if not sample_map[wes_sample] == plink_sample:
                    mismatches.append(l)
            # score >= 0.05 is not a good match - but check the top hit anyway just in case
            else:
                if sample_map[wes_sample] == plink_sample:
                    if plink_sample in plink_samples:
                        dodgy_samples[wes_sample] = {'line': l.rstrip(), 'status': 'match', 'in_plink': 'yes'}
                    else:
                        dodgy_samples[wes_sample] = {'line': l.rstrip(), 'status': 'match', 'in_plink': 'no'}
                else:
                    if plink_sample in plink_samples:
                        dodgy_samples[wes_sample] = {'line': l.rstrip(), 'status': 'mismatch', 'in_plink': 'yes'}
                    else:
                        dodgy_samples[wes_sample] = {'line': l.rstrip(), 'status': 'mismatch', 'in_plink': 'no'}

    with open(gtcheck_mismatches_file, 'w') as o1:
        o1.write("#samples with mismatches between WES ID and genotyping ID\n")
        o1.write(("\t").join(['#EGA', 'fam_id_top_hit', 'sum_discordance', 'nsites', "discordance/nsites"]))
        o1.write("\n")
        for m in mismatches:
            o1.write(m)

    with open(gtcheck_dodgy_samples_file, 'w') as o2:
        o2.write("#samples with high discordance score\n")
        o2.write(("\t").join(['#EGA', 'fam_id_top_hit', 'sum_discordance',
                 'nsites', "discordance/nsites", "status", "id_in_plink_file"]))
        o2.write("\n")
        for s in dodgy_samples.keys():
            o2.write(dodgy_samples[s]['line'] + "\t" + dodgy_samples[s]
                     ['status'] + "\t" + dodgy_samples[s]['in_plink'] + "\n")


def main():
    gtcheck_output_file = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/whole_exome_output/gtcheck_best_hits.txt'
    plinksamples_file = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/plink_samples.txt'
    metadata_file = '/lustre/scratch119/humgen/projects/birth_cohort_wes/genestack/all_samples_with_proceed_and_seq_info_and_warehouse_info.txt'
    gtcheck_dodgy_samples_file = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/whole_exome_output/gtcheck_samples_without_good_hit.txt'
    gtcheck_duplicates_file = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/whole_exome_output/gtcheck_duplicates.txt'
    gtcheck_mismatches_file = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/whole_exome_output/gtcheck_mismatches.txt'

    plink_samples = find_plink_samples(plinksamples_file)
    sample_map = parse_metadata(metadata_file, gtcheck_duplicates_file)
    parse_gtcheck_output(gtcheck_output_file, plink_samples, sample_map,
                         gtcheck_dodgy_samples_file, gtcheck_mismatches_file)


if __name__ == '__main__':
    main()
