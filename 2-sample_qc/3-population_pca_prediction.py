# population prediction with PCA
import hail as hl
import pyspark
import argparse
from gnomad.sample_qc.ancestry import assign_population_pcs
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
    

def create_1kg_mt(resourcedir: str, mtdir: str):
    '''
    Create matrixtable of 1kg data
    :param str resourcedir: resources directory
    :param str mtdir: matrixtable directory
    '''
    indir = resourcedir + "1kg_vcfs_filtered_by_wes_baits/"
    vcfheader = indir + "header.txt"
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


def run_pca(filtered_mt_file: str, pca_scores_file: str, pca_loadings_file: str, pca_evals_file: str):
    '''
    Run PCA before population prediction
    :param str filtered_mt_file: merged birth cohort wes and 1kg MT file annotated with pops and filtered
    :param str pca_sores_file: PCA scores HT file
    :param str pca_loadings_file: PCA scores HT file
    :param str pca_evals_file: PCA scores HT file
    '''
    mt = hl.read_matrix_table(filtered_mt_file)
    pca_evals, pca_scores, pca_loadings = hl.hwe_normalized_pca(mt.GT, k=10, compute_loadings=True)
    pca_scores = pca_scores.annotate(known_pop=mt.cols()[pca_scores.s].known_pop)
    pca_scores.write(pca_scores_file, overwrite=True)
    pca_loadings.write(pca_loadings_file, overwrite=True)
    with open(pca_evals_file, 'w') as f:
        for val in pca_evals:
            f.write(str(val) + "\n")



def predict_pops(pca_scores_file: str, pop_ht_file: str):
    '''
    Predict populations from PCA scores
    :param str pca_sores_file: PCA scores HT file
    :param str pop_ht_file: predicted population HT file
    '''
    pca_scores = hl.read_table(pca_scores_file)
    known_col = "known_pop"
    pop_ht, pop_clf = assign_population_pcs(pca_scores, pca_scores.scores, known_col=known_col, n_estimators=100, prop_train=0.8, min_prob=0.5)
    pop_ht.write(pop_ht_file, overwrite=True)
    # This is a work around for hail <0.2.88 - convert hail table to pandas df then run assign_population_pcs
    # pop_pca_scores = pca_scores.select(known_col, pca_scores=pca_scores.scores)
    # pop_pc_pd = pop_pca_scores.to_pandas()
    # pop_pc_pd = expand_pd_array_col(pop_pc_pd, "pca_scores", 10, 'PC')  

    # pop_pca_scores = pca_scores.select(known_col, pca_scores=pc_cols)
    # pop_pc_pd = pop_pca_scores.to_pandas()
    # pc_cols = [f"PC{i+1}" for i in range(10)]
    # pop_pd, pop_clf = assign_population_pcs(pop_pc_pd, pc_cols, known_col=known_col, n_estimators=100, prop_train=0.8, min_prob=0.5)


def main():
    #set up
    args = get_options()
    inputs = parse_config()
    mtdir = inputs['matrixtables_lustre_dir']
    resourcedir = inputs['resource_dir']
    mtdir2 = inputs['load_matrixtables_lustre_dir']

    #initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    #if needed, create 1kg matrixtable
    if args.kg_to_mt:
        create_1kg_mt(resourcedir, mtdir)

    #combine with 1KG data
    pruned_mt_file = mtdir + "mt_ldpruned.mt"
    merged_mt_file = mtdir + "merged_with_1kg.mt"
    if args.merge or args.run:
        merge_with_1kg(pruned_mt_file, mtdir, merged_mt_file)

    #annotate and filter
    filtered_mt_file = mtdir + "merged_with_1kg_filtered.mt"
    if args.filter or args.run:
        annotate_and_filter(merged_mt_file, resourcedir, filtered_mt_file)

    #run pca
    pca_scores_file = mtdir + "pca_scores_after_pruning.ht"
    pca_loadings_file = mtdir + "pca_loadings_after_pruning.ht"
    pca_evals_file = mtdir2 + "pca_evals_after_pruning.txt"#text file may need to be without file///
    if args.pca or args.run:
        run_pca(filtered_mt_file, pca_scores_file, pca_loadings_file, pca_evals_file)

    #assign pops
    pop_ht_file = mtdir + "pop_assignments.ht"
    if args.assign_pops or args.run:
        predict_pops(pca_scores_file, pop_ht_file)


if __name__ == '__main__':
    main()
