# 112 samples with high C>A resequenced to see if C>A artefact is fixed and to 
# see if this brings the nSNPs and rTiTv in line with what is expected

import hail as hl
import pyspark
from bokeh.plotting import output_file, save
from hail.plot import show, output_notebook
import pandas as pd
from wes_qc.utils.utils import parse_config


def load_vcfs_to_mt(indir, mtoutfile, tmp_dir, header):
    '''
    load VCFs and save as hail mt
    '''
    objects = hl.utils.hadoop_ls(indir)
    vcfs = [vcf["path"] for vcf in objects if (vcf["path"].startswith("file") and vcf["path"].endswith("vcf.gz"))]
    print("Loading VCFs")
    #create and save MT
    mt = hl.import_vcf(vcfs, array_elements_required=False, force_bgz=True, header_file = header)
    print("Saving as hail mt")
    mt.write(mtoutfile, overwrite=True)


def run_sample_qc_and_annotate(mtfile, spectrafile, plot_file):
    '''
    Run sample QC, annotate output with C>A fraction and plot nSNP vs rTiTv
    '''
    mt = hl.read_matrix_table(mtfile)
    mt =  hl.sample_qc(mt, name='sample_qc')
    ht = mt.cols()
    spec_df = pd.read_csv(spectrafile ,sep="\t")
    spec_df['CA'] = spec_df['C>A'] + spec_df['G>T']#C>A and G>T both count to the total fraction C>A
    spectra_ht = hl.Table.from_pandas(spec_df).key_by('sample')
    ht = ht.annotate(CA=spectra_ht[ht.s].CA)#annotate sample QC ht with CA fraction
    p = hl.plot.scatter(ht.sample_qc.n_snp, ht.sample_qc.r_ti_tv, xlabel='nSNP', ylabel='Ti/Tv ratio', label=ht.CA, collect_all = True)
    output_file(filename=plot_file)
    save(p)


def subset_mt_by_samples(mtinfile, sample_list, mtoutfile):
    '''
    Subset hail mt to just samples contained in a list
    '''
    mt = hl.read_matrix_table(mtinfile)
    sample_ht = hl.import_table(sample_list, no_header=True)
    sample_ht = sample_ht.key_by('f0')
    mt_filtered = mt.filter_cols(hl.is_defined(sample_ht[mt.s]))
    mt_filtered.write(mtoutfile, overwrite=True)


def main():
    #set up input variables
    inputs = parse_config()
    vcf_header = inputs['gatk_vcf_header']
    mtdir = inputs['load_matrixtables_lustre_dir']
    plot_dir = inputs['plots_dir_local']

    #initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    import_vcf_dir = "file:///lustre/scratch123/qc/repeat_samples_gatk/"

    #load VCFs into hail
    mtoutfile = mtdir + "gatk_unprocessed_repeat_samples.mt"
    load_vcfs_to_mt(import_vcf_dir, mtoutfile, tmp_dir, vcf_header)

    #run sample QC, annotated with fraction CA and plot nSNP/rTiTv
    spectrafile = "file:///lustre/scratch123/qc/repeat_samples_mutation_spectra/plots/proportions_per_person.txt"
    plot_file = plot_dir + "resequenced_samples_nSNP_rTiTv.html"
    run_sample_qc_and_annotate(mtoutfile, spectrafile, plot_file)

    #subset original mt to just the resequenced samples and repeat the run of sample qC, annotation with fraction CA in originals and plotting
    sample_list = "file:///lustre/scratch123/qc/mutation_spectra/samples.txt"
    mt_file_orig = mtdir +  "gatk_unprocessed.mt"
    mt_orig_subsetfile = mtdir + "resequenced_samples_original_data.mt"
    subset_mt_by_samples(mt_file_orig, sample_list, mt_orig_subsetfile)
    spectrafile_orig = "file:///lustre/scratch123/qc/repeat_samples_mutation_spectra/plots/proportions_per_person.txt"
    plot_file_orig = plot_dir + "resequenced_samples_nSNP_rTiTv_original_data.html"
    run_sample_qc_and_annotate(mt_orig_subsetfile, spectrafile_orig, plot_file_orig)

    #plot old vs new fraction C>A for each sample

if __name__ == '__main__':
    main() 