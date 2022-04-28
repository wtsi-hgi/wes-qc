# check the broad samples, using the Bristol metadata as the truth set
# output to contain:
# Broad ID
# Broad alternate ID (for those which are sometimes _2 )
# ALSPAC id
# Cram sent to Sanger
# EGA accession 1
# gtcheck_match_plink (does the top hit match the EGA accession in plink)
# gtcheck_match_wes (does the top hit match the EGA accession when WES data compared)
# EGA accession 2
# gtcheck_match_plink (does the top hit match the EGA accession in plink)
# gtcheck_match_wes (does the top hit match the EGA accession when WES data compared)

from nis import match


def parse_bristol_metadata(bristol_metadata_file):
    bristol_id_map = {}
    with open(bristol_metadata_file, 'r') as f:
        lines = f.readlines()
        for l in lines:
            ldata = l.split()
            if ldata[0] == '"Collaborator_Sample_ID"':
                continue
            broad_id = ldata[1].replace('"', '')
            alspac_id = ldata[0] + ldata[2].replace('"', '')
            bristol_id_map[broad_id] = {'alspac_id': alspac_id, 'in_plink': '',
                                        'alspac_plink_match': '', 'alternative_id': '', 'sent_to_sanger': 'no', 'ega_accs': {}}

    return bristol_id_map


def identify_sanger_only_samples(bristol_id_map, broad_wes_sample_file, additional_sanger_samples_file):
    '''
    identify samples which were sequences at sanger but not in the Bristol metadata file
    also identify alternate ids
    '''
    sanger_only = []
    with open(broad_wes_sample_file, 'r') as f:
        sequenced_samples = f.readlines()
        sequenced_samples = [x.rstrip() for x in sequenced_samples]

    for s in sequenced_samples:
        if not s in bristol_id_map.keys():
            s_mod = s[:-2]
            if s_mod in bristol_id_map.keys():
                bristol_id_map[s_mod]['alternative_id'] = s
            else:
                sanger_only.append(s)

    with open(additional_sanger_samples_file, 'w') as o:
        o.write(("\n").join(sanger_only))


def get_ega_accs(bristol_id_map, gtcheck_duplicates_file, ega_to_broad):
    '''
    Identify the broad samples that were also sequenced at Sanger, and their EGA accessions
    '''
    with open(gtcheck_duplicates_file, 'r') as f:
        lines = f.readlines()
        for l in lines:
            ldata = l.split()
            ega_to_broad[ldata[0]] = ldata[1]
            if ldata[1] in bristol_id_map.keys():
                bristol_id_map[ldata[1]]['ega_accs'][ldata[0]] = {'gtcheck_plink_match': '', 'gtcheck_wes_match': ''}
            elif ldata[1][:-2] in bristol_id_map.keys():  # also try alternate id
                altid = ldata[1][:-2]
                bristol_id_map[altid]['ega_accs'][ldata[0]] = {'gtcheck_plink_match': '', 'gtcheck_wes_match': ''}
            else:
                print(ldata[1] + " from duplicates file is not found in Bristol metadata")


def parse_wes_wes_gtcheck(bristol_id_map, gtcheck_wes_duplicates_file):
    '''
    Parse the output from the WES vs WES gtcheck for each broad sample
    '''
    with open(gtcheck_wes_duplicates_file, 'r') as f:
        lines = f.readlines()
        for l in lines:
            if l.startswith('#'):
                continue
            ldata = l.split('\t')
            ega = ldata[0]
            broad = ldata[5]
            status = ldata[6]
            if not ldata[7] == '':
                status = status + " " + ldata[7]
            if broad in bristol_id_map.keys():
                if ega in bristol_id_map[broad]['ega_accs'].keys():
                    bristol_id_map[broad]['ega_accs'][ega]['gtcheck_wes_match'] = status
                else:
                    print(ega + " not found in duplicates file")
            elif broad[:-2] in bristol_id_map.keys():  # check alternate id
                altid = broad[:-2]
                if ega in bristol_id_map[altid]['ega_accs'].keys():
                    bristol_id_map[altid]['ega_accs'][ega]['gtcheck_wes_match'] = status
                else:
                    print(ega + " not found in duplicates file")
            else:
                print(broad + " sample not found Bristol matadata")


def parse_plink_wes_mismatch(bristol_id_map, plink_gtcheck_mismatch_file, ega_to_broad):
    '''
    Parse the output from the WES ve plink gtcheck mismatches file for each sample. 
    If WES sample is sanger this needs to be mapped to the appropriate broad sample
    '''
    with open(plink_gtcheck_mismatch_file, 'r') as f:
        lines = f.readlines()
        for l in lines:
            if l.startswith('#'):
                continue
            ldata = l.split('\t')
            wes_id = ldata[0]
            if ldata[6] == 'yes':
                status = 'mismatch_expected_match_in_plink'
            elif ldata[6] == 'no':
                status = 'mismatch_expected_match_not_in_plink'
            else:
                status = ''
            if wes_id in bristol_id_map.keys():
                bristol_id_map[wes_id]['alspac_plink_match'] = status
            elif wes_id[:-2] in bristol_id_map.keys():
                altid = wes_id[:-2]
                bristol_id_map[altid]['alspac_plink_match'] = status
            elif wes_id in ega_to_broad.keys():
                broad_id = ega_to_broad[wes_id]
                if broad_id in bristol_id_map.keys():
                    if wes_id in bristol_id_map[broad_id]['ega_accs'].keys():
                        bristol_id_map[broad_id]['ega_accs'][wes_id]['gtcheck_plink_match'] = status
                    else:
                        print(wes_id + " in mismatches file not found in Bristol metadata/duplicates")
                elif broad_id[:-2] in bristol_id_map.keys():
                    altid = broad_id[:-2]
                    if wes_id in bristol_id_map[altid]['ega_accs'].keys():
                        bristol_id_map[altid]['ega_accs'][wes_id]['gtcheck_plink_match'] = status
                    else:
                        print(wes_id + " in mismatches file not found in Bristol metadata/duplicates")
                else:
                    print(broad_id + " mapping to " + wes_id + " in mismatches file not found in Bristol metadata")
            else:
                print(wes_id + " in mismatches file not found in Bristol metadata")


def parse_plink_wes_high_score(bristol_id_map, plink_gtcheck_high_score_file, ega_to_broad):
    '''
    Parse the output from the WES ve plink gtcheck high score file for each sample. 
    If WES sample is sanger this needs to be mapped to the appropriate broad sample
    '''
    with open(plink_gtcheck_high_score_file, 'r') as f:
        lines = f.readlines()
        for l in lines:
            if l.startswith('#'):
                continue
            ldata = l.split('\t')
            wes_id = ldata[0]
            if ldata[5] == 'match':
                status = 'match_high_score'
            else:
                if ldata[8] == 'yes':
                    status = 'mismatch_high_score'
                else:
                    status = 'mismatch_high_score_expected_match_not_in_plink'
            if wes_id in bristol_id_map.keys():
                bristol_id_map[wes_id]['alspac_plink_match'] = status
            elif wes_id[:-2] in bristol_id_map.keys():
                altid = wes_id[:-2]
                bristol_id_map[altid]['alspac_plink_match'] = status
            elif wes_id in ega_to_broad.keys():
                broad_id = ega_to_broad[wes_id]
                if broad_id in bristol_id_map.keys():
                    if wes_id in bristol_id_map[broad_id]['ega_accs'].keys():
                        bristol_id_map[broad_id]['ega_accs'][wes_id]['gtcheck_plink_match'] = status
                elif broad_id[:-2] in bristol_id_map.keys():
                    altid = broad_id[:-2]
                    if wes_id in bristol_id_map[altid]['ega_accs'].keys():
                        bristol_id_map[altid]['ega_accs'][wes_id]['gtcheck_plink_match'] = status


def parse_plink_matches(bristol_id_map, plink_gtcheck_match_file, ega_to_broad):
    '''
    Parse plink matches
    If WES sample is sanger this needs to be mapped to the appropriate broad sample
    '''
    with open(plink_gtcheck_match_file, 'r') as f:
        lines = f.readlines()
        for l in lines:
            if l.startswith('#'):
                continue
            ldata = l.split("\t")
            wes_id = ldata[0]
            status = 'match'
            if wes_id in bristol_id_map.keys():
                bristol_id_map[wes_id]['alspac_plink_match'] = status
            elif wes_id[:-2] in bristol_id_map.keys():
                altid = wes_id[:-2]
                bristol_id_map[altid]['alspac_plink_match'] = status
            elif wes_id in ega_to_broad.keys():
                broad_id = ega_to_broad[wes_id]
                if broad_id in bristol_id_map.keys():
                    if wes_id in bristol_id_map[broad_id]['ega_accs'].keys():
                        bristol_id_map[broad_id]['ega_accs'][wes_id]['gtcheck_plink_match'] = status
                elif broad_id[:-2] in bristol_id_map.keys():
                    altid = broad_id[:-2]
                    if wes_id in bristol_id_map[altid]['ega_accs'].keys():
                        bristol_id_map[altid]['ega_accs'][wes_id]['gtcheck_plink_match'] = status


def get_duplicate_info(bristol_id_map, broad_wes_sample_file, gtcheck_duplicates_file, gtcheck_wes_duplicates_file, plink_gtcheck_mismatch_file, plink_gtcheck_high_score_file, plink_gtcheck_match_file, additional_sanger_samples_file):
    # identify samples sequenced at sanger for which we have no metadata, and fill in alternate ids for when sanger and bristol id differ
    identify_sanger_only_samples(bristol_id_map, broad_wes_sample_file, additional_sanger_samples_file)
    ega_to_broad = {}
    get_ega_accs(bristol_id_map, gtcheck_duplicates_file, ega_to_broad)
    parse_wes_wes_gtcheck(bristol_id_map, gtcheck_wes_duplicates_file)
    parse_plink_wes_mismatch(bristol_id_map, plink_gtcheck_mismatch_file, ega_to_broad)
    parse_plink_wes_high_score(bristol_id_map, plink_gtcheck_high_score_file, ega_to_broad)
    parse_plink_matches(bristol_id_map, plink_gtcheck_match_file, ega_to_broad)


def print_outdata(bristol_id_map, outfile):
    '''
    Print outdata summary of all bristol samples
    '''
    with open(outfile, 'w') as o:
        header = ("\t").join(['Broad_id', 'Alternate_id', 'ALSPAC_id', 'in_plink', 'plink_gtcheck', 'sequenced_at_sanger', 'ega_acc_1',
                              'ega_acc_1_plink_gtcheck', 'ega_acc_1_wes_gtcheck', 'ega_acc_2', 'ega_acc_2_plink_gtcheck', 'ega_acc_2_wes_gtcheck'])
        o.write(header)
        o.write("\n")
        for s in bristol_id_map.keys():
            outdata = [s, bristol_id_map[s]['alternative_id'], bristol_id_map[s]['alspac_id'], bristol_id_map[s]['in_plink'], bristol_id_map[s]['sent_to_sanger']]
            if len(bristol_id_map[s]['ega_accs'].keys()) == 0:
                outdata = outdata + ['', '', '', '', '', '']
            elif len(bristol_id_map[s]['ega_accs'].keys()) == 1:
                for ega in bristol_id_map[s]['ega_accs'].keys():
                    outdata.append(ega)
                    outdata.append(bristol_id_map[s]['ega_accs'][ega]['gtcheck_plink_match'])
                    outdata.append(bristol_id_map[s]['ega_accs'][ega]['gtcheck_wes_match'])
                outdata = outdata + ['', '', '']

            elif len(bristol_id_map[s]['ega_accs'].keys()) == 2:
                    for ega in bristol_id_map[s]['ega_accs'].keys():
                        outdata.append(ega)
                        outdata.append(bristol_id_map[s]['ega_accs'][ega]['gtcheck_plink_match'])
                        outdata.append(bristol_id_map[s]['ega_accs'][ega]['gtcheck_wes_match'])
            else:
                print(s + " maps to >2 EGA accs, this shouldn't happen")
            o.write(("\t").join(outdata))
            o.write("\n")


def main():
    bristol_metadata_file = '/lustre/scratch119/humgen/projects/birth_cohort_wes/genestack/b3565_znumber.txt'
    broad_wes_sample_file = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/broad_sample_analysis/broad_wes_samples.txt'
    gtcheck_duplicates_file = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/check_wes_vs_wes/sample_pairs.txt'
    gtcheck_wes_duplicates_file = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/check_wes_vs_wes/gtcheck_analysis.txt'
    plink_gtcheck_mismatch_file = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/whole_exome_output/gtcheck_mismatches.txt'
    plink_gtcheck_high_score_file = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/whole_exome_output/gtcheck_samples_without_good_hit.txt'
    plink_gtcheck_match_file = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/whole_exome_output/gtcheck_matches.txt'
    # outfile for full info about each broad sample
    outfile = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/broad_sample_analysis/broad_samples_sanger_info.txt'
    # outfile for samples sequenced at snager which are not in the broad metadata
    additional_sanger_samples_file = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/broad_sample_analysis/broad_samples_sanger_not_in_metadata.txt'

    bristol_id_map = parse_bristol_metadata(bristol_metadata_file)
    get_duplicate_info(bristol_id_map, broad_wes_sample_file, gtcheck_duplicates_file, gtcheck_wes_duplicates_file,
                       plink_gtcheck_mismatch_file, plink_gtcheck_high_score_file, plink_gtcheck_match_file, additional_sanger_samples_file)
    print_outdata(bristol_id_map, outfile)


if __name__ == '__main__':
    main()
