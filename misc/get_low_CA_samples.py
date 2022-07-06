#filter the ALSPAC data to just samples with fraction CA < 0.08
#this is a stringent filter, the throey is that the high CA samples are creating excessive FP variants in variant QC
#using a very strict CA filter to test this

import hail as hl
import pyspark
import pandas as pd
from wes_qc.utils.utils import parse_config


def main():
    #set up input variables
    inputs = parse_config()
    mtdir = inputs['matrixtables_lustre_dir']

    spectrafile = "file:///lustre/scratch123/qc/mutation_spectra/plots/proportions_per_person.txt"
    #create pandas df from the spectra file, convert to hail table and filter
    spec_df = pd.read_csv(spectrafile ,sep="\t")
    spec_df['CA'] = spec_df['C>A'] + spec_df['G>T']#C>A and G>T both count to the total fraction C>A
    spectra_ht = hl.Table.from_pandas(spec_df).key_by('sample')
    spectra_ht = spectra_ht.filter(spectra_ht.CA < 0.08)#filter to only low CA samples
    #load in full MT
    mtfile = mtdir + "mt_pops_QC_filters_sequencing_location_and_superpop_sanger_only_after_sample_qc.mt"
    mt = hl.read_matrix_table(mtfile)
    #filter by sample
    filteredmt = mt.filter_cols(hl.is_defined(spectra_ht[mt.s]), keep=True)
    #save to file
    filteredmtfile = mtdir + "mt_varqc_splitmulti_lowCA_samples.mt"
    filteredmt.write(filteredmtfile, overwrite=True)


if __name__ == '__main__':
    main() 