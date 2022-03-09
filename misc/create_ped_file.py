# Create ped file for trios to use in variant QC
from wes_qc.utils.utils import parse_config


def parse_manifest(manifest_file):
    '''
    Parse manifest file and return a dict keyed by proband.
    We only want complete trios.
    '''
    trios = {}
    with open(manifest_file, 'r') as f:
        lines = f.readlines()
        for l in lines:
            if not l.startswith('sangersampleid'):
                linedata = l.split("\t")
                ega = linedata[25]
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
    manifest_file = resourcedir + "all_samples_with_proceed_and_seq_info_and_warehouse_info.txt"

    trios = triodata = parse_manifest(manifest_file)
    pedfile = resourcedir + "trios.ped"
    write_ped(trios, pedfile)

if __name__ == '__main__':
    main()