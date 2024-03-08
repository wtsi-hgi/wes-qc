import pandas
from io import StringIO
from typing import List, Union


class Family:
    parent_types = ('mother', 'father')

    def __init__(self, fid: str):
        self.id = fid
        self.mother = None
        self.father = None
        self.children = set()

    def set_parent(self, sample: str, parent_type: str):
        parent_type = parent_type.lower()
        assert parent_type in self.parent_types
        parent = getattr(self, parent_type)
        if parent is not None and parent != sample:
            raise ValueError(f'Conflict: self.{parent_type} = {parent} and {parent_type} = {sample}')
        setattr(self, parent_type, sample)

    def get_parent(self, parent_type: str) -> str:
        assert parent_type in self.parent_types
        parent = getattr(self, parent_type) or '0'
        return parent

    def add_child(self, sample: str):
        self.children.add(sample)

    def to_fam_records(self) -> List[List[str]]:
        records = []
        for child in self.children:
            record = [self.id, child, self.get_parent('father'), self.get_parent('mother'), '0', '-9']
            records.append(record)
        return records


def convert_dta(file: str):
    df = pandas.read_stata(file)
    buffer = StringIO()
    df.to_csv(buffer, sep=',', index=False)
    buffer.seek(0)
    return buffer


def parse_manifest(manifest_file: Union[str, StringIO]):
    trios = {}

    f = open(manifest_file) if isinstance(manifest_file, str) else manifest_file

    for line in f:
        sample1, sample2, relationship, type1, type2, pregnancy, preg_number = line.rstrip().split(',')

        if relationship == 'Partner' or relationship == 'relationship':
            continue

        if pregnancy in trios.keys():
            family = trios[pregnancy]
        else:
            family = Family(fid=pregnancy)
            trios[pregnancy] = family

        if relationship == 'Parent':
            family.set_parent(sample2, parent_type=type2)
            family.add_child(sample1)
        elif relationship == 'Child':
            family.set_parent(sample1, parent_type=type1)
            family.add_child(sample2)
        else:
            raise ValueError(f'Unknown relationship: {relationship}')

    return trios


def write_ped(trios, pedfile):
    with open(pedfile, 'w') as f:
        for trio in trios.values():
            for record in trio.to_fam_records():
                line = '\t'.join(record) + '\n'
                f.write(line)
        f.write('\n')


def replace_ids(input: str, idmap: str, output: str):
    df = pandas.read_table(input, sep='\t', header=None)
    mapping = pandas.read_csv(idmap, usecols=['BiBPersonID', 'SANGER SAMPLE ID'])
    mapping = mapping.append({'BiBPersonID': '0', 'SANGER SAMPLE ID': '0'}, ignore_index=True)

    missing_ids = []
    dt = df.merge(mapping, left_on=1, right_on='BiBPersonID', how='left').drop(columns='BiBPersonID').rename(columns={'SANGER SAMPLE ID': 'IID'})
    missing_ids.extend(dt.loc[pandas.isnull(dt.IID), 1])
    dt = dt.merge(mapping, left_on=2, right_on='BiBPersonID', how='left').drop(columns='BiBPersonID').rename(columns={'SANGER SAMPLE ID': 'FID'})
    missing_ids.extend(dt.loc[pandas.isnull(dt.FID), 2])
    dt = dt.merge(mapping, left_on=3, right_on='BiBPersonID', how='left').drop(columns='BiBPersonID').rename(columns={'SANGER SAMPLE ID': 'MID'})
    missing_ids.extend(dt.loc[pandas.isnull(dt.MID), 3])

    dt.fillna('0', inplace=True)
    dt = dt[[0, 'IID', 'FID', 'MID', 4, 5]]
    dt.to_csv(output, sep='\t', index=False, header=False)

    missing_ids_fn = output.replace('.ped', '.missing_ids.txt')
    missing_ids = list(set(missing_ids))
    with open(missing_ids_fn, 'w') as f:
        f.write('\n'.join(missing_ids))


def replace_another_ids(input: str, idmap: str, output: str):
    df = pandas.read_table(input, sep='\t', header=None)
    mapping = pandas.read_csv(idmap, sep='\t', header=None, names=['bibkidex', 'egan'])
    mapping = mapping.append({'bibkidex': '0', 'egan': '0'}, ignore_index=True)

    missing_ids = []
    dt = df.merge(mapping, left_on=1, right_on='bibkidex', how='left').drop(columns='bibkidex').rename(
        columns={'egan': 'IID'})
    missing_ids.extend(dt.loc[pandas.isnull(dt.IID), 1])
    dt = dt.merge(mapping, left_on=2, right_on='bibkidex', how='left').drop(columns='bibkidex').rename(
        columns={'egan': 'FID'})
    missing_ids.extend(dt.loc[pandas.isnull(dt.FID), 2])
    dt = dt.merge(mapping, left_on=3, right_on='bibkidex', how='left').drop(columns='bibkidex').rename(
        columns={'egan': 'MID'})
    missing_ids.extend(dt.loc[pandas.isnull(dt.MID), 3])

    dt.fillna('0', inplace=True)
    dt = dt[[0, 'IID', 'FID', 'MID', 4, 5]]
    dt.to_csv(output, sep='\t', index=False, header=False)

    missing_ids_fn = output.replace('.ped', '.missing_ids.txt')
    missing_ids = list(set(missing_ids))
    with open(missing_ids_fn, 'w') as f:
        f.write('\n'.join(missing_ids))


def main():
    manifest_file = '/lustre/scratch123/qc/BiB/Batch2_StudySampleIDs_RelatedPairs.csv'

    manifest_dta = '/lustre/scratch123/qc/BiB/BiB_CohortInfo__related_pairs.dta'
    manifest_file = convert_dta(manifest_dta)

    trios = parse_manifest(manifest_file)

    pedfile = "/lustre/scratch123/qc/BiB/trios.ped"
    write_ped(trios, pedfile)

    newpedfile = "/lustre/scratch123/qc/BiB/trios.bibkidex.ped"
    idmap = "/lustre/scratch123/qc/BiB/Batch2_StudySampleIDs_Matched.csv"
    replace_ids(pedfile, idmap, newpedfile)

    idmap = "/lustre/scratch123/qc/BiB/BiB.samples.tsv"
    anotherpedfile = "/lustre/scratch123/qc/BiB/trios.EGAN.ped"
    replace_another_ids(newpedfile, idmap, anotherpedfile)


if __name__ == '__main__':
    main()
