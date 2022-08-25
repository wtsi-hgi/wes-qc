# parse output of bcftools gtcheck
import gzip


def parse_gtcheck_file(gtcheck_file):
    '''
    parse gtcheck file and identify best hit
    '''
    best_hits = {}
    with open(gtcheck_file, 'r') as f:
        for line in f:
            if not line.startswith('DC'):
                continue
            linedata = line.split()
            sample = linedata[1]
            famid = linedata[2]
            discordance = float(linedata[3])
            nsites = int(linedata[5])
            if nsites > 0:
                score = discordance/nsites
            else:
                score = 1

            if sample in best_hits.keys():
                if score < best_hits[sample]['score']:
                    best_hits[sample] = {'fam': famid, 'discordance': discordance, 'nsites': nsites, 'score': score}
            else:
                best_hits[sample] = {'fam': famid, 'discordance': discordance, 'nsites': nsites, 'score': score}

    return best_hits


def write_output(outfile, best_hits):
    '''
    write output file
    '''
    print("Writing output")
    with open(outfile, 'w') as o:
        header = ("\t").join(['#EGA', 'broad_id_top_hit', 'sum_discordance', 'nsites', "discordance/nsites"])
        o.write(header + "\n")
        for sample in best_hits.keys():
            outline = ("\t").join([sample, best_hits[sample]['fam'], str(best_hits[sample]['discordance']), str(
                best_hits[sample]['nsites']), str(best_hits[sample]['score'])])
            o.write(outline + "\n")


def main():
    # gtcheck_output_file = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/whole_exome_output/gtcheck_top_50_matches.txt.gz'
    # outfile = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/whole_exome_output/gtcheck_best_hits.txt'
    gtcheck_output_file = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/compare_broad_sanger_vcfs/gtcheck_top_5_matches.txt'
    outfile = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/compare_broad_sanger_vcfs/gtcheck_best_hits.txt'
    best_hits = parse_gtcheck_file(gtcheck_output_file)
    write_output(outfile, best_hits)


if __name__ == '__main__':
    main()
