# parse output of bcftools gtcheck
import gzip
from os import listdir
from os.path import isfile, join


def read_output_files(gtcheck_output_dir):
    print('Reading GT check output files')
    output_data = {}
    gtcheck_files = [f for f in listdir(gtcheck_output_dir) if isfile(join(gtcheck_output_dir, f))]
    for gtf in gtcheck_files:
        gtf_path = gtcheck_output_dir + "/" + gtf
        print("Processing file: " + gtf_path)
        with gzip.open(gtf_path, 'rt') as zf:
            for line in zf:
                if line.startswith('DC'):
                    linedata = line.split()
                    ega = linedata[1]
                    famid = linedata[2]
                    discordance = float(linedata[3])
                    nsites = int(linedata[5])
                    idpair = ega + "_" + famid
                    if not idpair in output_data.keys():
                        output_data[idpair] = {'sum_discordance':0, 'nsites':0}
                    output_data[idpair]['sum_discordance'] +=  discordance
                    output_data[idpair]['nsites'] += nsites

    return output_data


def parse_output_data(output_data):
    '''
    identify best hit for each EGA/Z acc
    '''
    print("Identifying best matches for each sample")
    best_hits = {}
    for idpair in output_data.keys():
        idsplit = idpair.split("_")
        ega = idsplit[0]
        print(ega)
        fam = idsplit[1]
        dscore = output_data[idpair]['sum_discordance'] / output_data[idpair]['nsites']
        if ega not in best_hits.keys():
            best_hits[ega] = {'famid_best_hit': fam, 'sum_discordance_best_hit': output_data[idpair]['sum_discordance'],
                              'sum_nsites_best_hit': output_data[idpair]['nsites'], 'discordance_score_best_hit': dscore}
        else:
            if dscore < best_hits[ega]['discordance_score_best_hit']:
                #replace what is in best hits with current better scoring hit
                best_hits[ega]['sum_discordance_best_hit'] = output_data[idpair]['sum_discordance']
                best_hits[ega]['sum_nsites_best_hit'] = output_data[idpair]['nsites']
                best_hits[ega]['discordance_score_best_hit'] = dscore
                best_hits[ega]['famid_best_hit'] = fam


def write_output(best_hits, outfile):
    print("Writing output")
    with open(outfile, 'w') as o:
        header = ("\t").join(['EGA', 'fam_id_top_hit', 'sum_discordance', 'nsites', "discordance/nsites"])
        o.write(header)
        o.write("\n")
        for ega in best_hits.keys():
            outline = ("\t").join([ega, best_hits[ega]['famid_best_hit'], str(best_hits[ega]['sum_discordance_best_hit']), str(
                best_hits[ega]['sum_nsites_best_hit']), str(best_hits[ega]['discordance_score_best_hit'])])
            o.write(outline)
            o.write("\n")


def main():
    gtcheck_output_dir = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/per_chromosome_output/'
    output_data = read_output_files(gtcheck_output_dir)
    best_hits = parse_output_data(output_data)
    outfile = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/gtcheck_best_hits.txt'
    write_output(best_hits, outfile)

if __name__ == '__main__':
    main()
