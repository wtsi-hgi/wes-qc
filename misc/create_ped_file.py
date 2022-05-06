# Create ped file for trios to use in variant QC
import gzip
from wes_qc.utils.utils import parse_config



def get_samples_to_exclude(sample_qc_fails, gtcheck_mismatches):
    '''
    Get samples to exclude from trios file - those which fail sample QC or those where 
    bcftools gtcheck has identified a mismatch
    '''
    to_exclude = {}
    with gzip.open(sample_qc_fails,'r') as q:
        for line in q:
            l = line.decode("utf-8")
            if l.startswith('EGAN'):
                l = l.rstrip()
                to_exclude[l] = 1

    with open(gtcheck_mismatches, 'r') as m:
        lines = m.readlines()
        for l in lines:
            if l.startswith('EGAN'):
                ldata = l.split()
                to_exclude[ldata[0]] = 1

    return to_exclude
        


def parse_manifest(manifest_file, samples_to_exclude):
    '''
    Parse manifest file and return a dict keyed by proband.
    We only want complete trios.
    '''
    trios = {}
    count = 0
    with open(manifest_file, 'r') as f:
        lines = f.readlines()
        for l in lines:
            if not l.startswith('sangersampleid'):
                linedata = l.split("\t")
                ega = linedata[25]
                if ega in samples_to_exclude.keys():
                    count += 1
                    print('excl ' + str(count))
                    continue
                persontype = linedata[4]
                famid = linedata[8]
                trio = linedata[6]
                sex = linedata[26]
                sexcode = '1'
                if sex == 'Female':
                    sexcode = '2'
                elif sex == 'Male':
                    sexcode = '1'
                else:# in ALSPAC there is one proband with unknown gender in manifest, removing this
                    continue

                if not ega.startswith('EGA'):#get rid of any people without a valid EGA acc
                    continue
                if trio == 'yes':
                    # fix A or B
                    persontype = persontype.replace(' ','')
                    if famid not in trios.keys():
                        trios[famid] = {}
                    trios[famid][persontype] = {}
                    trios[famid][persontype]['ega'] = ega
                    trios[famid][persontype]['sex'] = sex

    return trios


def write_ped(trios, pedfile):
    peddata = []
    for famid in trios.keys():
        if 'M' in trios[famid].keys() and 'P' in trios[famid].keys():#both parents present
            for person in trios[famid].keys():
                if person not in ['M', 'P']:
                    pedline = ("\t").join([ famid, trios[famid][person]['ega'], trios[famid]['P']['ega'], trios[famid]['M']['ega'], trios[famid][person]['sex'], '0' ])
                    peddata.append(pedline)

    with open(pedfile, 'w') as o:
        o.write(("\n").join(peddata))


def main():
    inputs = parse_config()
    resourcedir = inputs['resource_dir_local']
    annotdir = inputs['annotation_lustre_dir_local']
    manifest_file = resourcedir + "all_samples_with_proceed_and_seq_info_and_warehouse_info.txt"
    sample_qc_fails = annotdir + "samples_failing_qc.tsv.bgz"
    gtcheck_mismatches = resourcedir + "gtcheck_mismatches.txt"

    samples_to_exclude = get_samples_to_exclude(sample_qc_fails, gtcheck_mismatches)
    trios  = parse_manifest(manifest_file, samples_to_exclude)
    pedfile = resourcedir + "trios.ped"
    write_ped(trios, pedfile)

if __name__ == '__main__':
    main()