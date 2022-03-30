# parse output of bcftools gtcheck


from email.quoprimime import header_check
from random import betavariate
from tkinter import E


def read_output_files(output_file_list):
    output_data = {}
    with open(output_file_list, 'r') as f:
        lines = f.readlines()
        for l in lines:
            gtfile = l.rstrip()
            with open(gtfile, 'r') as g:
                glines = g.readlines()
                for gl in glines:
                    if gl.startswith('DC'):
                        gldata = gl.split()
                        ega = gldata[1]
                        famid = gldata[2]
                        discordance = float(gldata[3])
                        nsites = int(gldata[5])
                        if not ega in output_data.keys():
                            output_data[ega] = {}
                        if not famid in output_data[ega].keys():
                            output_data[ega][famid] = {'sum_discordance':0, 'nsites':0}
                        output_data[ega][famid]['sum_discordance'] +=  discordance
                        output_data[ega][famid]['nsites'] += nsites

    return output_data


def parse_output_data(output_data):
    best_hits = {}
    for ega in output_data.keys():
        best_discordance_score = 100.0
        best_sum_discordance = 1.0
        best_sum_nsites = 0
        best_hit = 'NA'
        for famid in output_data[ega].keys():
            dscore = output_data[ega][famid]['sum_discordance'] / output_data[ega][famid]['nsites']
            if dscore < best_discordance_score:
                best_discordance_score = dscore
                best_sum_discordance = output_data[ega][famid]['sum_discordance']
                best_sum_nsites = output_data[ega][famid]['nsites']
                best_hit = famid
        best_hits[ega] = {'famid_best_hit':best_hit, 'sum_discordance_best_hit': best_sum_discordance, 'sum_nsites_best_hit': best_sum_nsites, 'discordance_score_best_hit': best_discordance_score }

    return best_hits

def write_output(best_hits, outfile):
    with open(outfile, 'w') as o:
        header = ("\t").join(['EGA', 'fam_id_top_hit', 'sum_discordance', 'nsites', "discordance/nsites"])
        o.write(header)
        o.write("\n")
        for ega in best_hits.keys():
            outline = ("\t").join([ ega, best_hits[ega]['famid_best_hit'], str(best_hits[ega]['sum_discordance_best_hit']), str(best_hits[ega]['sum_nsites_best_hit']), str(best_hits[ega]['discordance_score_best_hit']) ])
            o.write(outline)
            o.write("\n")


def main():
    output_file_list = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/gtcheck_output_files.txt'
    output_data = read_output_files(output_file_list)
    best_hits = parse_output_data(output_data)
    outfile = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/gtcheck_best_hits.txt'
    write_output(best_hits, outfile)

if __name__ == '__main__':
    main() 