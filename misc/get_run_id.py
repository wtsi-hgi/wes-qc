#get sequencing run IDs for Sanger samples to annotate plots for sample QC
#run on hgi-farm5 not on os tenant due to scratch119 files

sequencing_info_file = '/lustre/scratch119/humgen/projects/birth_cohort_wes/genestack/6407_QC_pass_seq_Batches.csv'
idmap_file = '/lustre/scratch119/humgen/projects/birth_cohort_wes/genestack/sequencescape_dump_study_6407.txt'
outfile = '/lustre/scratch123/hgi/projects/birth_cohort_wes/qc/resources/sequencing_batches.txt'

sanger_ega = {}
ega_runid = {}

with open(idmap_file, 'r') as f:
    lines = f.readlines()
    for l in lines:
        if not l.startswith('sanger_sample_id'):
            ldata = l.split()
            sanger_ega[ldata[0]] = ldata[2]

with open(sequencing_info_file, 'r') as f2:
    lines2 = f2.readlines()
    for l2 in lines2:
        if not l2.startswith('sanger_sample_id'):
            l2data = l2.split(",")
            ega = sanger_ega[l2data[0]]
            runid = l2data[3]
            ega_runid[ega] = runid

with open(outfile, 'w') as o:
    o.write(("\t").join(['ega', 'runid']))
    o.write("\n")
    for e in ega_runid.keys():
        o.write(("\t").join([e, ega_runid[e]]))
        o.write("\n")
