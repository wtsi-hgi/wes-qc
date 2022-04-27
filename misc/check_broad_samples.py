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
            print(ldata)
            exit(0)




def main():
    bristol_metadata_file = '/lustre/scratch119/humgen/projects/birth_cohort_wes/genestack/b3565_znumber.txt'
    broad_sample_file = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/broad_sample_analysis/broad_wes_samples.txt'

    bristol_id_map = parse_bristol_metadata(bristol_metadata_file)



if __name__ == '__main__':
    main()