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

def parse_bristol_metadata(bristol_metadata_file):
    bristol_id_map = {}
    with open(bristol_metadata_file, 'r') as f:
        lines = f.readlines()
        for l in lines:
            ldata = l.split()
            if ldata[0] == '"Collaborator_Sample_ID"':
                continue
            broad_id = ldata[1].replace('"','')
            alspac_id = ldata[0] + ldata[2].replace('"','')
            bristol_id_map[broad_id] = alspac_id
    
    return bristol_id_map


def identify_sanger_only_samples(bristol_id_map, broad_wes_sample_file, additional_sanger_samples_file):
    sanger_only = []
    with open(broad_wes_sample_file, 'r') as f:
        sequenced_samples = f.readlines()
        sequenced_samples = [x.rstrip() for x in sequenced_samples]

    for s in sequenced_samples:
        if not s in bristol_id_map.keys():
            sanger_only.append(s)
    
    with open(additional_sanger_samples_file, 'w') as o:
        o.write(("/n").join(sanger_only))


def get_duplicate_info(bristol_id_map, broad_wes_sample_file, gtcheck_duplicates_file, gtcheck_wes_duplicates_file, outfile, additional_sanger_samples_file):
    #identify samples sequenced at sanger for which we have no metadata
    identify_sanger_only_samples(bristol_id_map, broad_wes_sample_file, additional_sanger_samples_file)


def main():
    bristol_metadata_file = '/lustre/scratch119/humgen/projects/birth_cohort_wes/genestack/b3565_znumber.txt'
    broad_wes_sample_file = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/broad_sample_analysis/broad_wes_samples.txt'
    gtcheck_duplicates_file = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/whole_exome_output/gtcheck_duplicates.txt'
    gtcheck_wes_duplicates_file = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/check_wes_vs_wes/gtcheck_analysis.txt'
    outfile = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/broad_sample_analysis/broad_samples_sanger_info.txt'#outfile for full info about each broad sample
    additional_sanger_samples_file = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/broad_sample_analysis/broad_samples_sanger_not_in_metadata.txt'#outfile for samples sequenced at snager which are not in the broad metadata

    bristol_id_map = parse_bristol_metadata(bristol_metadata_file)
    get_duplicate_info(bristol_id_map, broad_wes_sample_file, gtcheck_duplicates_file, gtcheck_wes_duplicates_file, outfile, additional_sanger_samples_file)



if __name__ == '__main__':
    main()