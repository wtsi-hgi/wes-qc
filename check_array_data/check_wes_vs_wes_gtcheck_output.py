# take file with the best hit for each sample from gtcheck
# check that the identifier matches that in sample_apirs.txt and that the score is good (< 0.05)
# write to an output file

def get_sample_pairs(expected_pairs_file):
    pairs = {}
    with open(expected_pairs_file, 'r') as f:
        lines = f.readlines()
        for l in lines:
            ldata = l.split()
            pairs[ldata[0]] = ldata[1]
    return pairs


def analyse_gtcheck_output(gtcheck_output_file, expected_pairs, output_file):
    outdata = []
    with open(gtcheck_output_file, 'r') as f:
        lines = f.readlines()
        for l in lines:
            if l.startswith('#EGA'):
                outline = l.rstrip() + "\t" + "expected_match" + "\t" + "match_status" + "\t" + "note" + "\n"
                outdata.append(outline)
            else:
                l = l.rstrip()
                ldata = l.split()
                match = "no"
                expected_match = expected_pairs[ldata[0]]
                if ldata[1] == expected_match:
                    match = 'yes'
                note = ''
                if float(ldata[4]) > 0.05:
                    note = 'high_score'
                outline = l + "\t" + expected_match + "\t" + match + "\t" + note + "\n"
                outdata.append(outline)
        
    with open(output_file, 'w') as o:
        for l in outdata:
            o.write(l)


def main():
    gtcheck_output_file = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/check_wes_vs_wes/gtcheck_best_hits.txt'
    expected_pairs_file = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/check_wes_vs_wes/sample_pairs.txt'
    output_file = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/check_wes_vs_wes/gtcheck_analysis.txt'

    expected_pairs = get_sample_pairs(expected_pairs_file)
    analyse_gtcheck_output(gtcheck_output_file, expected_pairs, output_file)


if __name__ == '__main__':
    main()