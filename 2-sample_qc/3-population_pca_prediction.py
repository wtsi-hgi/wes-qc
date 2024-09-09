# population prediction with PCA
import hail as hl
import pyspark
import os
import argparse
import pandas as pd
from typing import Optional
from gnomad.sample_qc.ancestry import assign_population_pcs
from utils.utils import parse_config
from utils.config import path_local, path_spark

def get_options():
    '''
    Get options from the command line
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("-k", "--kg_to_mt", 
        help="convert 1kg data to matrixtable", action="store_true")
    parser.add_argument("-m", "--merge",
        help="merge alspac mt with 1kg mt", action="store_true")
    parser.add_argument("-f", "--filter",
        help="annotate and filter merged mt", action="store_true")
    parser.add_argument("-p", "--pca",
        help="run pca", action="store_true")
    parser.add_argument("-a", "--assign_pops",
        help="assign populations", action="store_true")
    parser.add_argument("-r", "--run",
        help="run all steps except kg_to_mt", action="store_true")

    args = parser.parse_args()

    return args
    

def create_1kg_mt(config: dict):
    '''
    Create matrixtable of 1kg data
    :param str resourcedir: resources directory
    :param str mtdir: matrixtable directory
    '''
    conf = config['step2']['create_1kg_mt']

    # TODO: correct resource dir: /lustre/scratch126/WES_QC/resources/1000g_VCFs
    indir = path_spark(conf['indir'])
    vcfheader = conf['vcfheader'] # TODO: DEBUG: not used during tests on small dataset
    objects = hl.utils.hadoop_ls(indir)
    vcfs = [vcf["path"] for vcf in objects if (vcf["path"].startswith("file") and vcf["path"].endswith("vcf.gz"))]
    print("Loading VCFs")
    #create and save MT
    # TODO: make header optional (don't have one for mini 1000g test)
    mt = hl.import_vcf(vcfs, array_elements_required=False, force_bgz=True)
    print("Saving as hail mt")
    mt_out_file = path_spark(conf['mt_out_file'])
    mt.write(mt_out_file, overwrite=True)

    return mt


def merge_with_1kg(pruned_mt_file: hl.MatrixTable, config: dict) -> hl.MatrixTable:
    '''
    merge the birth cohort WES ld-pruned data with 1kg data of known population
    :param str pruned_mt_file: ld pruned MT file
    :param str mtdir: resources directory
    :param str merged_mt_file: merged output MT file
    '''
    conf = config['step2']['merge_with_1_kg']
    print("Merging with 1kg data")
    # mt = hl.read_matrix_table(pruned_mt_file)
    mt = pruned_mt_file # FIXME: just testing accepting mt instead of a path
    kg_mt_file = path_spark(conf['kg_mt_file'])
    kg_mt = hl.read_matrix_table(kg_mt_file)
    # in order to create a union dataset the entries and column fields must be 
    # the same in each dataset. The following 2 lines take care of this.
    kg_mt = kg_mt.select_entries(kg_mt.GT)
    mt = mt.drop('callrate', 'f_stat', 'is_female') # TODO: move to config or not? seems to be specific to the qc pipeline so no need to make variable
    # union cols gives all samples and the rows which are found in both
    mt_joined = mt.union_cols(kg_mt)
    merged_mt_file = path_spark(conf['merged_mt_outfile'])
    mt_joined.write(merged_mt_file, overwrite=True)

    return mt_joined


def annotate_and_filter(merged_mt_file: hl.MatrixTable, config: dict) -> hl.MatrixTable:
    '''
    Annotate with known pops for 1kg samples and filter to remove long range LD regions, 
    rare variants, palidromic variants, low call rate and HWE filtering
    :param str merged_mt_file: merged birth cohort wes and 1kg MT file
    :param str resourcedir: resources directory
    :param str filtered_mt_file: merged birth cohort wes and 1kg MT file annotated with pops and filtered
    '''
    conf = config['step2']['annotate_and_filter']
    print("Adding population annotation for 1kg samples")
    # mt = hl.read_matrix_table(merged_mt_file)
    mt = merged_mt_file

    # The following annotates by population
    # pops_file = resourcedir + "integrated_call_samples.20130502.ALL.ped"
    # cohorts_pop = hl.import_table(pops_file, delimiter="\t").key_by('Individual ID')
    # mt = mt.annotate_cols(known_pop=cohorts_pop[mt.s].Population)

    # The following is 1kg superpop

    # TODO: remove leading slash
    pops_file = path_spark(conf['pops_file'])
    cohorts_pop = hl.import_table(pops_file, delimiter="\t").key_by('Sample name')
    mt = mt.annotate_cols(known_pop=cohorts_pop[mt.s]['Superpopulation code'])

    print("Filtering variants")
    mt_vqc = hl.variant_qc(mt, name='variant_QC_Hail')
    mt_vqc_filtered = mt_vqc.filter_rows(
        (mt_vqc.variant_QC_Hail.call_rate >= conf['call_rate']) &
        (mt_vqc.variant_QC_Hail.AF[1] >= conf['AF']) &
        (mt_vqc.variant_QC_Hail.p_value_hwe >= conf['p_value_hwe'])
    )
    # TODO: this file must be .bed (?)
    long_range_ld_file = path_spark(conf['long_range_ld_file'])
    long_range_ld_to_exclude = hl.import_bed(long_range_ld_file, reference_genome='GRCh38')
    mt_vqc_filtered = mt_vqc_filtered.filter_rows(hl.is_defined(long_range_ld_to_exclude[mt_vqc_filtered.locus]), keep=False)
    mt_non_pal = mt_vqc_filtered.filter_rows((mt_vqc_filtered.alleles[0] == "G") & (mt_vqc_filtered.alleles[1] == "C"), keep=False)
    mt_non_pal = mt_non_pal.filter_rows((mt_non_pal.alleles[0] == "C") & (mt_non_pal.alleles[1] == "G"), keep=False)
    mt_non_pal = mt_non_pal.filter_rows((mt_non_pal.alleles[0] == "A") & (mt_non_pal.alleles[1] == "T"), keep=False)
    mt_non_pal = mt_non_pal.filter_rows((mt_non_pal.alleles[0] == "T") & (mt_non_pal.alleles[1] == "A"), keep=False)

    filtered_mt_file = path_spark(conf['filtered_mt_outfile'])
    mt_non_pal.write(filtered_mt_file, overwrite=True)

    return mt_non_pal


def run_pca(filtered_mt_file: hl.MatrixTable, config: dict):
    '''
    Run PCA before population prediction
    :param str filtered_mt_file: merged birth cohort wes and 1kg MT file annotated with pops and filtered
    :param str pca_sores_file: PCA scores HT file
    :param str pca_loadings_file: PCA scores HT file
    :param str pca_evals_file: PCA scores HT file
    '''
    conf = config['step2']['run_pca']
    # mt = hl.read_matrix_table(filtered_mt_file)
    mt = filtered_mt_file # FIXME: testing accepting mt
    #exclude EGAN00004311029 this is a huge outlier on PCA and skews all the PCs
    to_exclude = ['EGAN00004311029'] # TODO: add as a list parameter to the config
    mt = mt.filter_cols(~hl.set(to_exclude).contains(mt.s)) 

    # TODO: move k to config
    pca_evals, pca_scores, pca_loadings = hl.hwe_normalized_pca(mt.GT, k=4, compute_loadings=True)
    pca_scores = pca_scores.annotate(known_pop=mt.cols()[pca_scores.s].known_pop)
    pca_scores_file = path_spark(conf['pca_scores_outfile'])
    pca_scores.write(pca_scores_file, overwrite=True)
    pca_loadings_file = path_spark(conf['pca_loadings_outfile'])
    pca_loadings.write(pca_loadings_file, overwrite=True)
    pca_evals_file = path_local(conf['pca_evals_outfile'])

    with open(pca_evals_file, 'w') as f:
        for val in pca_evals:
            f.write(str(val) + "\n")

    # don't return anything as there are multiple outputs


def append_row(df, row):
    return pd.concat([
                df, 
                pd.DataFrame([row], columns=row.index)]
           ).reset_index(drop=True)


def predict_pops(config: dict):
    '''
    Predict populations from PCA scores
    :param str pca_sores_file: PCA scores HT file
    :param str pop_ht_file: predicted population HT file
    :param str pop_ht_tsv: population tsv file
    '''
    conf = config['step2']['predict_pops']
    pca_scores_file = path_spark(conf['pca_scores_file'])
    pca_scores = hl.read_table(pca_scores_file)
    known_col = "known_pop"
    pop_ht, pop_clf = assign_population_pcs(pca_scores, pca_scores.scores, 
                                            known_col=known_col, 
                                            n_estimators=conf['gnomad_pc_n_estimators'], 
                                            prop_train=conf['gnomad_prop_train'], 
                                            min_prob=conf['gnomad_min_prob'])
    pop_ht_file = path_spark(conf['pop_ht_outfile'])
    pop_ht.write(pop_ht_file, overwrite=True)
    #convert to pandas and put in only pops files, add excluded sample back
    pop_ht_df = pop_ht.to_pandas()
    pop_ht_df2 =pop_ht_df[['s', 'pop']]
    new_row = pd.Series({'s':'EGAN00004311029', 'pop':'oth'})
    pop_ht_df2 = append_row(pop_ht_df2, new_row)
    print(pop_ht_df2)

    pop_ht_tsv = path_local(conf['pop_ht_outtsv'])
    print(pop_ht_tsv)

    pop_ht_df2.to_csv(pop_ht_tsv, sep = '\t')
    
    # This is a work around for hail <0.2.88 - convert hail table to pandas df then run assign_population_pcs
    # pop_pca_scores = pca_scores.select(known_col, pca_scores=pca_scores.scores)
    # pop_pc_pd = pop_pca_scores.to_pandas()
    # pop_pc_pd = expand_pd_array_col(pop_pc_pd, "pca_scores", 10, 'PC')  

    # pop_pca_scores = pca_scores.select(known_col, pca_scores=pc_cols)
    # pop_pc_pd = pop_pca_scores.to_pandas()
    # pc_cols = [f"PC{i+1}" for i in range(10)]
    # pop_pd, pop_clf = assign_population_pcs(pop_pc_pd, pc_cols, known_col=known_col, n_estimators=100, prop_train=0.8, min_prob=0.5)


def main():
    # get cli args
    args = get_options()
    config = parse_config()
    tmp_dir = config['general']['tmp_dir']
    mtdir = config['general']['matrixtables_dir']

    #initialise hail
    # tmp_dir = "hdfs://spark-master:9820/"
    # sc = pyspark.SparkContext()
    sc = pyspark.SparkContext.getOrCreate()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38", idempotent=True)

    #if needed, create 1kg matrixtable
    if args.kg_to_mt:
        # the output is not used in this main func
        created_1kg_mt = create_1kg_mt(config)

    #combine with 1KG data
    if args.merge or args.run:
        pruned_mt_file = path_spark(os.path.join(mtdir, 'mt_ldpruned.mt'))
        pruned_mt = hl.read_matrix_table(pruned_mt_file)

        # the output is not used in this main func
        mt_joined = merge_with_1kg(pruned_mt, config)

    #annotate and filter
    if args.filter or args.run:
        merged_mt_file = path_spark(config['step2']['merge_with_1_kg']['merged_mt_outfile'])
        merged_mt = hl.read_matrix_table(merged_mt_file)
        
        # the output is not used in this main func
        mt_non_pal = annotate_and_filter(merged_mt, config)

    #run pca
    if args.pca or args.run:
        filtered_mt_file = path_spark(config['step2']['annotate_and_filter']['filtered_mt_outfile'])
        filtered_mt = hl.read_matrix_table(filtered_mt_file)
        run_pca(filtered_mt, config)

    #assign pops
    if args.assign_pops or args.run:
        predict_pops(config)

if __name__ == '__main__':
    main() 
