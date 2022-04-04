#get sequencing run IDs for Sanger samples to annotate plots for sample QC

sequencing_info_file = '/lustre/scratch119/humgen/projects/birth_cohort_wes/genestack/6407_QC_pass_seq_Batches.csv'
idmap_file = '/lustre/scratch119/humgen/projects/birth_cohort_wes/genestack/sequencescape_dump_study_6407.txt'

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
        print(l2)
        exit(0)