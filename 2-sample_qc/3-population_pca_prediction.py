# population prediction with PCA
import hail as hl
import pyspark
import argparse
import pandas as pd
from gnomad.sample_qc.ancestry import assign_population_pcs, pc_project
from utils.utils import parse_config

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
    

def create_1kg_mt(resourcedir: str, mtdir: str, indir: str, vcfheader: str):
    '''
    Create matrixtable of 1kg data
    :param str resourcedir: resources directory
    :param str mtdir: matrixtable directory
    '''
    objects = hl.utils.hadoop_ls(indir)
    vcfs = [vcf["path"] for vcf in objects if (vcf["path"].startswith("file") and vcf["path"].endswith("vcf.gz"))]
    print("Loading VCFs")
    #create and save MT
    mt = hl.import_vcf(vcfs, array_elements_required=False, force_bgz=True, header_file = vcfheader)
    print("Saving as hail mt")
    mt_out_file = mtdir + "kg_wes_regions.mt"
    mt.write(mt_out_file, overwrite=True)


def merge_with_1kg(pruned_mt_file: str, mtdir: str, merged_mt_file: str):
    '''
    merge the birth cohort WES ld-pruned data with 1kg data of known population
    :param str pruned_mt_file: ld pruned MT file
    :param str mtdir: resources directory
    :param str merged_mt_file: merged output MT file
    '''
    print("Merging with 1kg data")
    mt = hl.read_matrix_table(pruned_mt_file)
    kg_mt_file = mtdir + "kg_wes_regions.mt"
    kg_mt = hl.read_matrix_table(kg_mt_file)
    # in order to create a union dataset the entries and column fields must be 
    # the same in each dataset. The following 2 lines take care of this.
    kg_mt = kg_mt.select_entries(kg_mt.GT)
    mt = mt.drop('callrate', 'f_stat', 'is_female')
    # union cols gives all samples and the rows which are found in both
    mt_joined = mt.union_cols(kg_mt)
    mt_joined.write(merged_mt_file, overwrite=True)


def annotate_and_filter(merged_mt_file: str, resourcedir: str, filtered_mt_file: str):
    '''
    Annotate with known pops for 1kg samples and filter to remove long range LD regions, 
    rare variants, palidromic variants, low call rate and HWE filtering
    :param str merged_mt_file: merged birth cohort wes and 1kg MT file
    :param str resourcedir: resources directory
    :param str filtered_mt_file: merged birth cohort wes and 1kg MT file annotated with pops and filtered
    '''
    print("Adding population annotation for 1kg samples")
    mt = hl.read_matrix_table(merged_mt_file)   

    # The following annotates by population
    # pops_file = resourcedir + "integrated_call_samples.20130502.ALL.ped"
    # cohorts_pop = hl.import_table(pops_file, delimiter="\t").key_by('Individual ID')
    # mt = mt.annotate_cols(known_pop=cohorts_pop[mt.s].Population)

    # The following is 1kg superpop
    pops_file = resourcedir + "/igsr_samples.tsv"
    cohorts_pop = hl.import_table(pops_file, delimiter="\t").key_by('Sample name')
    mt = mt.annotate_cols(known_pop=cohorts_pop[mt.s]['Superpopulation code'])

    print("Filtering variants")
    mt_vqc = hl.variant_qc(mt, name='variant_QC_Hail')
    mt_vqc_filtered = mt_vqc.filter_rows(
        (mt_vqc.variant_QC_Hail.call_rate >= 0.99) &
        (mt_vqc.variant_QC_Hail.AF[1] >= 0.05) &
        (mt_vqc.variant_QC_Hail.p_value_hwe >= 1e-5)
    )
    long_range_ld_file = resourcedir + "long_range_ld_regions_chr.txt"
    long_range_ld_to_exclude = hl.import_bed(long_range_ld_file, reference_genome='GRCh38')
    mt_vqc_filtered = mt_vqc_filtered.filter_rows(hl.is_defined(long_range_ld_to_exclude[mt_vqc_filtered.locus]), keep=False)
    mt_non_pal = mt_vqc_filtered.filter_rows((mt_vqc_filtered.alleles[0] == "G") & (mt_vqc_filtered.alleles[1] == "C"), keep=False)
    mt_non_pal = mt_non_pal.filter_rows((mt_non_pal.alleles[0] == "C") & (mt_non_pal.alleles[1] == "G"), keep=False)
    mt_non_pal = mt_non_pal.filter_rows((mt_non_pal.alleles[0] == "A") & (mt_non_pal.alleles[1] == "T"), keep=False)
    mt_non_pal = mt_non_pal.filter_rows((mt_non_pal.alleles[0] == "T") & (mt_non_pal.alleles[1] == "A"), keep=False)

    mt_non_pal.write(filtered_mt_file, overwrite=True)


def run_pca(filtered_mt_file: str, pca_scores_file: str, pca_1kg_scores_file: str, pca_1kg_loadings_file: str, pca_1kg_evals_file: str):
    '''
    Run PCA before population prediction
    :param str filtered_mt_file: merged birth cohort wes and 1kg MT file annotated with pops and filtered
    :param str pca_sores_file: PCA scores HT file
    :param str pca_loadings_file: PCA scores HT file
    :param str pca_evals_file: PCA scores HT file
    '''
    mt = hl.read_matrix_table(filtered_mt_file)
    #divide matrix to make a projection
    mt_kg = mt.filter_cols(hl.is_defined(mt.known_pop))
    mt_study = mt.filter_cols(hl.is_missing(mt.known_pop))
    #PCA for 1000 Genomes
    pca_evals, pca_scores, pca_loadings = hl.hwe_normalized_pca(mt_kg.GT, k=10, compute_loadings=True)
    pca_scores = pca_scores.annotate(known_pop=mt_kg.cols()[pca_scores.s].known_pop)
    pca_af_ht = mt_kg.annotate_rows(pca_af=hl.agg.mean(mt_kg.GT.n_alt_alleles()) / 2).rows()
    pca_loadings = pca_loadings.annotate(pca_af=pca_af_ht[pca_loadings.key].pca_af)
    #saving files
    pca_scores.write(pca_1kg_scores_file, overwrite=True)
    pca_loadings.write(pca_1kg_loadings_file, overwrite=True)
    with open(pca_1kg_evals_file, 'w') as f:
        for val in pca_evals:
            f.write(str(val) + "\n")
    pca_scores = pca_scores.drop(pca_scores.known_pop)
    #projection of samples on precomputed PCs and combining of two PCA_scores tables
    projection_PCA_scores=pc_project(mt_study, pca_loadings, loading_location='loadings', af_location='pca_af')
    union_PCA_scores=pca_scores.union(projection_PCA_scores)
    union_PCA_scores = union_PCA_scores.annotate(known_pop=mt.cols()[pca_scores.s].known_pop)
    #saving union_PCA_scores
    union_PCA_scores.write(pca_scores_file, overwrite=True)


def append_row(df, row):
    return pd.concat([
                df, 
                pd.DataFrame([row], columns=row.index)]
           ).reset_index(drop=True)


def predict_pops(pca_scores_file: str, pop_ht_file: str, pop_ht_tsv):
    '''
    Predict populations from PCA scores
    :param str pca_sores_file: PCA scores HT file
    :param str pop_ht_file: predicted population HT file
    :param str pop_ht_tsv: population tsv file
    '''
    pca_scores = hl.read_table(pca_scores_file)
    known_col = "known_pop"
    pop_ht, pop_clf = assign_population_pcs(pca_scores, pca_scores.scores, known_col=known_col, n_estimators=100, prop_train=0.8, min_prob=0.5)
    pop_ht.write(pop_ht_file, overwrite=True)
    #convert to pandas and put in only pops files, add excluded sample back
    pop_ht_df = pop_ht.to_pandas()
    pop_ht_df2 =pop_ht_df[['s', 'pop']]
    #new_row = pd.Series({'s':'EGAN00004311029', 'pop':'oth'})
    #pop_ht_df2 = append_row(pop_ht_df2, new_row)
    print(pop_ht_df2)
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

def filter_matrix (mt: hl.MatrixTable, long_range_ld_file: str):
    #use only autosomes
    mt=mt.filter_rows(mt.locus.in_autosome())
    #split multiallelic variants and remove them
    mt = hl.split_multi_hts(mt)#this shouldn't do anything as only biallelic sites are used
    mt = mt.filter_rows(mt.was_split==True, keep=False)
    #keep only SNPs
    mt = mt.filter_rows(hl.is_snp(mt.alleles[0], mt.alleles[1]))
    #keep good variants using hail variant_qc and thre filters
    mt_vqc = hl.variant_qc(mt, name='variant_QC_Hail')
    mt_vqc_filtered = mt_vqc.filter_rows(
        (mt_vqc.variant_QC_Hail.call_rate >= 0.99) &
        (mt_vqc.variant_QC_Hail.AF[1] >= 0.05) &
        (mt_vqc.variant_QC_Hail.p_value_hwe >= 1e-5)
    )
    #remove long ld regions
    long_range_ld_to_exclude = hl.import_bed(long_range_ld_file, reference_genome='GRCh38')
    mt_vqc_filtered = mt_vqc_filtered.filter_rows(hl.is_defined(long_range_ld_to_exclude[mt_vqc_filtered.locus]), keep=False)
    #remove palindromes
    mt_non_pal = mt_vqc_filtered.filter_rows((mt_vqc_filtered.alleles[0] == "G") & (mt_vqc_filtered.alleles[1] == "C"), keep=False)
    mt_non_pal = mt_non_pal.filter_rows((mt_non_pal.alleles[0] == "C") & (mt_non_pal.alleles[1] == "G"), keep=False)
    mt_non_pal = mt_non_pal.filter_rows((mt_non_pal.alleles[0] == "A") & (mt_non_pal.alleles[1] == "T"), keep=False)
    mt_non_pal = mt_non_pal.filter_rows((mt_non_pal.alleles[0] == "T") & (mt_non_pal.alleles[1] == "A"), keep=False)
    
    return mt_non_pal

def merge_and_prune (mt_file: str, kg_mt_file: str, long_range_ld_file: str, pops_file: str, merged_mt_file: str, pruned_mt_file: str):
    #filter both matrices
    mt=hl.read_matrix_table(mt_file)
    kg_mt=hl.read_matrix_table(kg_mt_file)
    mt_filtered=filter_matrix(mt, long_range_ld_file)
    kg_filtered=filter_matrix(kg_mt, long_range_ld_file)
    #removing and adding needed entries to replicate filtered_mt_file structure 
    mt_filtered=mt_filtered.drop('AD', 'DP', 'GQ', 'MIN_DP', 'PGT', 'PID', 'PL', 'PS', 'SB', 'RGQ', 'callrate')
    kg_filtered = kg_filtered.select_entries(kg_filtered.GT)
    #merging matrices
    mt_merged = mt_filtered.union_cols(kg_filtered)
    #saving the merged matrix
    mt_merged=mt_merged.checkpoint(merged_mt_file, overwrite=True)
    #prunning
    pruned_ht = hl.ld_prune(mt_merged.GT, r2=0.2)
    pruned_mt = mt_merged.filter_rows(hl.is_defined(pruned_ht[mt_merged.row_key]))
    pruned_mt = pruned_mt.select_entries(GT=hl.unphased_diploid_gt_index_call(pruned_mt.GT.n_alt_alleles()))
    #annotating known poppulations
    cohorts_pop = hl.import_table(pops_file, delimiter="\t").key_by('Sample name')
    pruned_mt = pruned_mt.annotate_cols(known_pop=cohorts_pop[pruned_mt.s]['Superpopulation code'])
    #saving matrix
    pruned_mt.write(pruned_mt_file, overwrite=True)

def main():
    #set up
    args = get_options()
    inputs = parse_config()
    mtdir = inputs['matrixtables_lustre_dir']
    resourcedir = inputs['resource_dir']
    mtdir2 = inputs['load_matrixtables_lustre_dir']

    #initialise hail
    #tmp_dir = "hdfs://spark-master:9820/"
    tmp_dir = inputs['tmp_dir']
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    #if needed, create 1kg matrixtable
    if args.kg_to_mt:
        kgdir = inputs['kg_dir']
        vcfheader = inputs['kg_vcfheader']
        create_1kg_mt(resourcedir, mtdir, kgdir, vcfheader)


    #combine with 1KG data
    #pruned_mt_file = mtdir + "mt_ldpruned.mt"
    mt_file = mtdir + "mt_hard_filters_annotated.mt"
    mt_kg_file = mtdir + "kg_wes_regions.mt"
    merged_mt_file = mtdir + "merged_with_1kg.mt"
    merged_pruned_mt_file = mtdir + "merged_with_1kg.pruned.mt"
    long_range_ld_file = resourcedir + "long_ld_regions.hg38.bed"
    pops_file = resourcedir + "igsr_samples.tsv"
    if args.merge or args.run:
        merge_and_prune (mt_file, mt_kg_file, long_range_ld_file, pops_file, merged_mt_file, merged_pruned_mt_file)

    #annotate and filter
    #filtered_mt_file = mtdir + "merged_with_1kg_filtered.mt"
    #if args.filter or args.run:
    #    annotate_and_filter(merged_mt_file, resourcedir, filtered_mt_file)

    #run pca
    pca_scores_file = mtdir + "pca_scores_after_pruning.test.ht"
    pca_1kg_scores_file = mtdir + "pca_scores_1kg.test.ht"
    pca_1kg_loadings_file = mtdir + "pca_loadings_1kg.test.ht"
    pca_1kg_evals_file = mtdir2 + "pca_evals_1kg.test.txt"#text file may need to be without file///

    if args.pca or args.run:
        run_pca(merged_pruned_mt_file, pca_scores_file, pca_1kg_scores_file, pca_1kg_loadings_file, pca_1kg_evals_file)

    #assign pops
    pop_ht_file = mtdir + "pop_assignments.ht"
    pop_ht_tsv = mtdir2 + "pop_assignemtnts.tsv"
    if args.assign_pops or args.run:
        predict_pops(pca_scores_file, pop_ht_file, pop_ht_tsv)


if __name__ == '__main__':
    main() 
