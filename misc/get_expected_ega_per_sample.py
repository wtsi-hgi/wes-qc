def get_z_to_a(infofile):
    z_to_a = {}
    with open(infofile, 'r') as f:
        lines = f.readlines()
        for l in lines:
            if l.startswith('Broad_id'):
                continue
            ldata = l.split('\t')
            zid = ldata[0]
            altzid = ldata[1]
            aid = ldata[2]
            z_to_a[zid] = aid
            if not altzid == '':
                z_to_a[altzid] = aid
    return z_to_a


def get_a_to_ega(infofile):
    a_to_e = {}
    with open(infofile, 'r') as f:
        lines = f.readlines()
        for l in lines:
            if l.startswith('Broad_id'):
                continue
            ldata = l.split('\t')
            aid = ldata[2]
            ega = ldata[5]
            a_to_e[aid] = ega
    return a_to_e


def get_ids_with_multiple_egas(infofile):
    ids = {}
    with open(infofile, 'r') as f:
        lines = f.readlines()
        for l in lines:
            if l.startswith('Broad_id'):
                continue
            ldata = l.split('\t')
            ega2 = ldata[8]
            aid = ldata[2]
            if not ega2 == '':
                ids[aid] = 1
    return ids


def get_problem_familes(infofile):
    problem_fams = {}
    with open(infofile, 'r') as f:
        lines = f.readlines()
        for l in lines:
            if l.startswith('Broad_id'):
                continue
            ldata = l.split('\t')
            aid = ldata[2]
            if ldata[5] == '':
                continue
            if not ldata[6] == 'match' or not ldata[7] == 'yes':
                problem_fams[aid] = 1
    return problem_fams


def main():
    match_info_file = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/compare_broad_sanger_vcfs/broad_samples_sanger_info.txt'
    samples_file = '/lustre/scratch119/humgen/projects/birth_cohort_wes/alspac/existing_data/WES_from_Broad/samples.txt'
    z_to_a_ids = get_z_to_a(match_info_file)#map ids in form Z1234 to 1234A, include alternates eg Z1234_2 
    a_to_ega_ids = get_a_to_ega(match_info_file)#map ida from 1234A to EGA accession, some 1234A idsa have no ega accession - these will be blank
    aid_multiple_egas = get_ids_with_multiple_egas(match_info_file)#get family ids (1234A) of individuals with >1 EGA - just print these out for now and check manually
    #flag up family ids with potential problems - column ega_acc_1_wes_gtcheck is not equal to yes (column 
    problem_families = get_problem_familes(match_info_file)

    outfile = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/compare_broad_sanger_vcfs/sample_mapping,txt'
    with open(outfile, 'w') as o:
        with open(samples_file, 'r') as f:
            lines = f.readlines()
            for l in lines:
                mapline = ''
                sid = l.rstrip()
                if sid == '':
                    continue
                #print(sid)
                if sid in a_to_ega_ids.keys():
                    mapline = sid + "\t" + a_to_ega_ids[sid] + "\n"
                    if sid in aid_multiple_egas.keys():
                        print(sid + " associated with multiple EGAs")
                    if sid in problem_families.keys():
                        print(sid + " may have id match problems")
                elif sid in z_to_a_ids.keys():
                    sid2 = z_to_a_ids[sid]
                    mapline = sid + "\t" + a_to_ega_ids[sid2] + "\n"
                    if sid2 in aid_multiple_egas.keys():
                        print(sid2 + " associated with multiple EGAs")
                    if sid2 in problem_families.keys():
                        print(sid2 + " may have id match problems")
                else:
                    maplines = sid + "\t" + "NA" + "\n"
                o.write(mapline)


if __name__ == '__main__':
    main()
