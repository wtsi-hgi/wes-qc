# run bcftools gtcheck on all samples in parallel over all shards
import argparse

def get_options():
    '''
    Get options from the command line
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--run", 
        help="run  gtcheck", action="store_true")    
    parser.add_argument("-p", "--parse",
        help="concatenate and parse output", action="store_true")
    args = parser.parse_args()

    return args


def create_ids_file(wdir, metadata_file):
    '''
    create sample id mapping file
    '''
    outdata = []
    with open(metadata_file, 'r') as m:
        lines = m.readlines()
        for l in linedata:
            linedata = lines.split("\t")
            print(linedata)
            exit(0)


def submit_gtcheck_jobs():
    '''
    find all gtcheck vcfs
    submit gtcheck jobs to farm - one per vcf
    '''
    pass


def concatenate_outputs():
    '''
    concatenate outputs, add up total discordance and total sites for each sample
    identifer pair and produce file with ega_id, discordance, sites and discorance/sites
    '''
    pass


def main():
    wdir = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/run_getcheck/'
    metadata_file = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/resources/all_samples_with_proceed_and_seq_info_and_warehouse_info_egas_from_mlwh.txt'

    args = get_options()

    if args.run:
        create_ids_file(wdir, metadata_file)
        submit_gtcheck_jobs()

    if args.parse:
        concatenate_outputs()

    if not args.parse and not args.run:
        print("Either --run (-r) or --parse (-p) is required")


if __name__ == '__main__':
    main()