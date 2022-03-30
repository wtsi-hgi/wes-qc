# parse output of bcftools gtcheck


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
                        print(gldata)
                        exit(0)


    return output_data


def main():
    output_file_list = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/check_array_genotypes/gtcheck_output_files.txt'
    output_data = read_output_files(output_file_list)

if __name__ == '__main__':
    main() 