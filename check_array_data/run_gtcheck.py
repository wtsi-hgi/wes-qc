# run bcftools gtcheck on all samples in parallel over all shards
import argparse
import subprocess

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


def runcommand(cmd):
    try:
        byteoutput = subprocess.check_output(cmd, shell=True)
        return byteoutput.decode('UTF-8').rstrip()
    except subprocess.CalledProcessError as e:
        print(e.output)
        return "Error in command"


def create_ids_file(wdir, metadata_file, gatk_vcf_dir, plink_vcf):
    '''
    create sample id mapping file
    check that the EGA accession is in the EGA VCFs (if it is in one it will be in all)
    check that the family ID is in the plink VCFs
    if one or other missing then skip that pair and report
    '''
    outdata = {}
    with open(metadata_file, 'r') as m:
        lines = m.readlines()
        for l in lines:
            linedata = l.split("\t")
            if linedata[0] == 'sangersampleid':
                continue
            ega = linedata[25]
            familyid = linedata[11]
            if not familyid.endswith('A or B'):
                outdata[ega] = familyid
    
    #run bcftools query -l to get list of plink vcf samples
    bcftools_cmd_p = "bcftools query -l " + plink_vcf
    plink_samples_str = runcommand(bcftools_cmd_p)
    plink_vcf_samples = plink_samples_str.split()
    print(plink_vcf_samples)
    exit(0)


    gatk_vcf_samples = []
    #run bcftools query -l to get list of gatk vcf samples
    gatk_vcf = gatk_vcf_dir + "chr10_101587143_chr10_101849993.gatk.vcf.gz"
    bcftools_cmd_g = "bcftools query -l " + gatk_vcf


    # for each sample pair - check if it is in both VCFs, if so add to the output 
    # file, if not report to command line




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
    gatk_vcf_dir = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/gatk_vcfs/'
    plink_vcf = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/plink_vcf_b38_liftover.vcf.gz'

    args = get_options()

    if args.run:
        create_ids_file(wdir, metadata_file, gatk_vcf_dir, plink_vcf)
        submit_gtcheck_jobs()

    if args.parse:
        concatenate_outputs()

    if not args.parse and not args.run:
        print("Either --run (-r) or --parse (-p) is required")


if __name__ == '__main__':
    main()