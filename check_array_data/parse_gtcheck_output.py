# parse output of bcftools gtcheck
import gzip
import os
import subprocess
from os import listdir
from os.path import isfile, join

def get_sample(sample_list, jobindex):
    '''
    get sample corresponding to the job index
    '''
    sample_num = jobindex - 1
    with open(sample_list, 'r') as f:
        lines = f.readlines()
        sample = lines[sample_num]
        sample = sample.rstrip()

    return sample


def runcommand(cmd):
    try:
        byteoutput = subprocess.check_output(cmd, shell=True)
        return byteoutput.decode('UTF-8').rstrip()
    except subprocess.CalledProcessError as e:
        print(e.output)
        return "Error in command"
        

def read_gtcheck_outputs(sample, gtcheck_output_dir):
    '''
    Extract data from one sample from the output files and combine
    '''
    print('Extracting information from GT check output files for ' + sample)
    output_data = {}
    gtcheck_files = [f for f in listdir(gtcheck_output_dir) if isfile(join(gtcheck_output_dir, f))]
    for gtf in gtcheck_files:
        gtf_path = gtcheck_output_dir + "/" + gtf
        print("Processing file: " + gtf_path)
        grepcmd = "zgrep " + sample + " " + gtf_path
        grepout = runcommand(grepcmd)
        greplines = grepout.split('\n')
        for gl in greplines:
            gldata = gl.split()
            ega = gldata[1]
            famid = gldata[2]
            discordance = float(gldata[3])
            nsites = int(gldata[5])
            idpair = ega + "_" + famid
            if not idpair in output_data.keys():
                output_data[idpair] = {'sum_discordance':0, 'nsites':0}
            output_data[idpair]['sum_discordance'] +=  discordance
            output_data[idpair]['nsites'] += nsites

    return output_data


# def read_output_files(gtcheck_output_dir):
#     print('Reading GT check output files')
#     output_data = {}
#     gtcheck_files = [f for f in listdir(gtcheck_output_dir) if isfile(join(gtcheck_output_dir, f))]
#     for gtf in gtcheck_files:
#         gtf_path = gtcheck_output_dir + "/" + gtf
#         print("Processing file: " + gtf_path)
#         with gzip.open(gtf_path, 'rt') as zf:
#             for line in zf:
#                 if line.startswith('DC'):
#                     linedata = line.split()
#                     ega = linedata[1]
#                     famid = linedata[2]
#                     discordance = float(linedata[3])
#                     nsites = int(linedata[5])
#                     idpair = ega + "_" + famid
#                     if not idpair in output_data.keys():
#                         output_data[idpair] = {'sum_discordance':0, 'nsites':0}
#                     output_data[idpair]['sum_discordance'] +=  discordance
#                     output_data[idpair]['nsites'] += nsites

#     return output_data



# def parse_output_data(output_data):
#     '''
#     identify best hit for each EGA/Z acc
#     '''
#     print("Identifying best matches for each sample")
#     best_hits = {}
#     for idpair in output_data.keys():
#         idsplit = idpair.split("_")
#         ega = idsplit[0]
#         fam = idsplit[1]
#         dscore = output_data[idpair]['sum_discordance'] / output_data[idpair]['nsites']
#         if ega not in best_hits.keys():
#             best_hits[ega] = {'famid_best_hit': fam, 'sum_discordance_best_hit': output_data[idpair]['sum_discordance'],
#                               'sum_nsites_best_hit': output_data[idpair]['nsites'], 'discordance_score_best_hit': dscore}
#         else:
#             if dscore < best_hits[ega]['discordance_score_best_hit']:
#                 #replace what is in best hits with current better scoring hit
#                 best_hits[ega]['sum_discordance_best_hit'] = output_data[idpair]['sum_discordance']
#                 best_hits[ega]['sum_nsites_best_hit'] = output_data[idpair]['nsites']
#                 best_hits[ega]['discordance_score_best_hit'] = dscore
#                 best_hits[ega]['famid_best_hit'] = fam

def parse_output_data(output_data):
    '''
    identify best hit for each sample
    '''
    print("Identifying best matches for each sample")
    best_hit = {'fam_id':'NA', 'dscore':'NA', 'sum_discordance':'NA', 'sum_nsites':'NA'}
    for idpair in output_data.keys():
        idsplit = idpair.split("_")
        fam = idsplit[1]
        dscore = output_data[idpair]['sum_discordance'] / output_data[idpair]['nsites']
        if best_hit['fam_id'] == 'NA':
            best_hit['fam_id'] = fam
            best_hit['dscore'] = dscore
            best_hit['sum_discordance'] = output_data[idpair]['sum_discordance']
            best_hit['sum_nsites'] = output_data[idpair]['nsites']
        elif dscore < best_hit['dscore']:
            best_hit['fam_id'] = fam
            best_hit['dscore'] = dscore
            best_hit['sum_discordance'] = output_data[idpair]['sum_discordance']
            best_hit['sum_nsites'] = output_data[idpair]['nsites']
    
    return best_hit


def write_header(outfile):
     with open(outfile, 'w') as o:
        header = ("\t").join(['EGA', 'fam_id_top_hit', 'sum_discordance', 'nsites', "discordance/nsites"])
        o.write(header)


def write_output(sample, best_hit, outfile):
    print("Writing output")
    with open(outfile, 'a') as o:
        outline = ("\t").join([sample, best_hit['fam_id'], str(best_hit['sum_discordance']), str(best_hit['sum_nsites']), str(best_hit['dscore']) ])
        o.write(outline + "\n")

# def write_output(best_hits, outfile):
#     print("Writing output")
#     with open(outfile, 'w') as o:
#         header = ("\t").join(['EGA', 'fam_id_top_hit', 'sum_discordance', 'nsites', "discordance/nsites"])
#         o.write(header)
#         o.write("\n")
#         for ega in best_hits.keys():
#             outline = ("\t").join([ega, best_hits[ega]['famid_best_hit'], str(best_hits[ega]['sum_discordance_best_hit']), str(
#                 best_hits[ega]['sum_nsites_best_hit']), str(best_hits[ega]['discordance_score_best_hit'])])
#             o.write(outline)
#             o.write("\n")


def main():
    gtcheck_output_dir = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/per_chromosome_output/'
    sample_list = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/gatk_samples.txt'
    outfile = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/gtcheck_best_hits.txt'
    jobindex = int(os.environ['LSB_JOBINDEX'])
    if jobindex == 1:#only write the header for the first job index
        write_header(outfile)

    sample = get_sample(sample_list, jobindex)
    #output_data = read_output_files(gtcheck_output_dir)
    output_data = read_gtcheck_outputs(sample, gtcheck_output_dir)
    best_hit = parse_output_data(output_data)
    write_output(sample, best_hit, outfile)

if __name__ == '__main__':
    main()
