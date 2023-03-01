# Create ped file for trios to use in variant QC
import gzip
import os.path
import logging
import pandas
from utils.utils import parse_config


def get_trios_to_exclude(unrelated_parents_file):
    '''
    Get trios to exclude due to unrelated parents
    '''
    to_exclude = {}
    with open(unrelated_parents_file, 'r') as f:
        lines = f.readlines()
        for l in lines:
            if not l.startswith('family'):
                ldata = l.split()
                to_exclude[ldata[0]] = 1

    return to_exclude


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
        

def parse_manifest(manifest_file, samples_to_exclude, trios_to_exclude):
    '''
    Parse manifest file and return a dict keyed by proband.
    We only want complete trios and we exclude samples which fail QC
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
                    continue
                persontype = linedata[4]
                famid = linedata[8]
                trio = linedata[6]
                sex = linedata[26]
                sexcode = '1'
                trioid = famid + persontype
                if trioid in trios_to_exclude.keys():
                    continue
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
                    trios[famid][persontype]['sex'] = sexcode

    return trios


def write_ped(trios, pedfile):
    peddata = []
    for famid in trios.keys():
        if 'M' in trios[famid].keys() and 'P' in trios[famid].keys():#both parents present
            for person in trios[famid].keys():
                if person not in ['M', 'P']:
                    fam = famid + person#this allows for families with >1 child
                    pedline = ("\t").join([ fam, trios[famid][person]['ega'], trios[famid]['P']['ega'], trios[famid]['M']['ega'], trios[famid][person]['sex'], '0' ])
                    peddata.append(pedline)

    with open(pedfile, 'w') as o:
        o.write(("\n").join(peddata))


class FamilyError(Exception):
    def __init__(self, df: pandas.DataFrame):
        self.df = df
        fam_id = df.hurles_fid.unique()[0]
        message = f'Bad family detected: {fam_id}'
        logging.error(message)
        super(FamilyError, self).__init__(message)


def metadata_to_ped(filename: str) -> pandas.DataFrame:
    metadata = pandas.read_table(filename, sep='\t')
    ped = []
    failed_data = []
    for family_id, family in metadata.groupby('hurles_fid'):
        try:
            row = family_to_ped(family, fid=family_id)
            ped.append(row)
        except FamilyError as e:
            failed_data.append(e.df)

    ped = pandas.concat(ped)
    return ped


def family_to_ped(df: pandas.DataFrame, fid: str) -> pandas.DataFrame:
    ped = []

    mother = get_person(df, person_type='M')
    father = get_person(df, person_type='F')
    if len(mother) + len(father) != 2:
        raise FamilyError(df)

    children = get_person(df, person_type='C')
    for child in children:
        if child == '0':
            continue
        gender = df.loc[df.sanger_sample_id == child, 'sex'].values
        ped_row = [fid, child, father[0], mother[0], gender[0], '-9']
        ped.append(ped_row)

    ped = pandas.DataFrame.from_records(ped)
    return ped


def get_person(df: pandas.DataFrame, person_type: str):
    assert person_type in ('M', 'F', 'C')
    if person_type not in df.mfc.values:
        return ['0']
    person_id = df.loc[df.mfc == person_type, 'sanger_sample_id']
    person_gender = df.loc[df.sanger_sample_id.isin(person_id), 'sex']
    if person_type == 'M' and any(person_gender != 2):
        logging.error(f'Wrong gender {person_gender} for mother {person_id}')
        raise FamilyError(df)
    if person_type == 'F' and any(person_gender != 1):
        logging.error(f'Wrong gender {person_gender} for father {person_id}')
        raise FamilyError(df)
    return person_id.values


def exclude_failed_samples(ped: pandas.DataFrame, failed_qc_file: str) -> pandas.DataFrame:
    df = pandas.read_csv(failed_qc_file, compression='gzip')
    ped = ped[~ped[1].isin(df.s)]
    ped[2].replace(df.s.to_list(), value='0', inplace=True)
    ped[3].replace(df.s.to_list(), value='0', inplace=True)
    return ped


def main():
    inputs = parse_config()
    resourcedir = inputs['resource_dir_local']
    annotdir = inputs['annotation_lustre_dir_local']
    manifest_file = resourcedir + "all_samples_with_proceed_and_seq_info_and_warehouse_info.txt"
    sample_qc_fails = annotdir + "samples_failing_qc.tsv.bgz"
    gtcheck_mismatches = resourcedir + "gtcheck_mismatches.txt"
    unrelated_parents_file = annotdir + "unrelated_parents.txt"
    metadata_file = os.path.join(resourcedir, 'GDAC_2021_03_HURLES_Sanger_mcs_basic_demographics_v0003_shareable_20220215.txt')
    pedfile = annotdir + "trios.ped"

    ped = metadata_to_ped(metadata_file)
    ped.to_csv(annotdir + "trios_all.ped", header=False, sep='\t', index=False)

    ped = exclude_failed_samples(ped, sample_qc_fails)
    ped.to_csv(pedfile, header=False, sep='\t', index=False)

    # trios_to_exclude = get_trios_to_exclude(unrelated_parents_file)
    # samples_to_exclude = get_samples_to_exclude(sample_qc_fails, gtcheck_mismatches)
    # trios = parse_manifest(manifest_file, samples_to_exclude, trios_to_exclude)
    # write_ped(trios, pedfile)


if __name__ == '__main__':
    main()
