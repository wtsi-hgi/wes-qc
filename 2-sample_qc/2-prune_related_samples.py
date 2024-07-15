# identify and prune related samples prior to PCA
# use mt with hard filters and sex annotation from 2-sample_qc/1-hard_filters_sex_annotation.py
import hail as hl
import pyspark
from wes_qc.utils.utils import parse_config
from bokeh.plotting import save, output_file


def prune_mt(mtin: hl.MatrixTable, mtoutfile: str):
    '''
    Splits multiallelic sites and runs ld pruning
    Filter to autosomes before LD pruning to decrease sample size - autosomes only wanted for later steps
    :param MatrixTable mtin: input MT containing variants to be pruned
    :param str mtdir: directory output matrix tables are written to
    '''
    print("Filtering to autosomes")
    mtin = mtin.filter_rows(mtin.locus.in_autosome())
    print("Splitting multiallelic sites")
    # mtin = hl.split_multi(mtin)
    mtin = hl.split_multi_hts(mtin)#this shouldn't do anything as only biallelic sites are used
    print("Performing LD pruning")
    # TODO: r2 to config?
    pruned_ht = hl.ld_prune(mtin.GT, r2=0.2)
    pruned_mt = mtin.filter_rows(hl.is_defined(pruned_ht[mtin.row_key]))
    pruned_mt = pruned_mt.select_entries(GT=hl.unphased_diploid_gt_index_call(pruned_mt.GT.n_alt_alleles()))
    pruned_mt.write(mtoutfile, overwrite=True)


def run_pc_relate(pruned_mt_file: str, relatedness_ht_file: str, samples_to_remove_file: str, scores_file: str):
    '''
    Runs PC relate on pruned MT
    :param str pruned_mt_file: ld pruned MT file
    :param str relatedness_ht_file: relatedness ht file
    :param str samples_to_remove_file: file samples to remove ht is written to
    :param str scores_file: file scores ht is written to
    '''
    print("Running PC relate")
    pruned_mt = hl.read_matrix_table(pruned_mt_file)
    # TODO: DEBUG: k=10 didn't work for the small test dataset, throwing an error:
    # "Error summary: IllegalArgumentException: requirement failed: Dimension mismatch!: _a.cols == _b.rows (11 != 10)"
    eig, scores, _ = hl.hwe_normalized_pca(pruned_mt.GT, k=3, compute_loadings=False)
    scores.write(scores_file, overwrite=True)

    print("Calculating relatedness")
    relatedness_ht = hl.pc_relate(pruned_mt.GT, min_individual_maf=0.05, scores_expr=scores[pruned_mt.col_key].scores, block_size=4096, min_kinship=0.05, statistics='kin2')
    relatedness_ht.write(relatedness_ht_file, overwrite=True)
    #prune individuals to be left with unrelated - creates a table containing one column - samples to remove
    pairs = relatedness_ht.filter(relatedness_ht['kin'] > 0.125)
    related_samples_to_remove = hl.maximal_independent_set(pairs.i, pairs.j, keep=False)
    related_samples_to_remove.write(samples_to_remove_file, overwrite=True)   


def run_population_pca(pruned_mt_file: str, pca_mt_file: str, mtdir: str, plotdir: str, samples_to_remove_file: str):
    '''
    Runs PCA and creates a matrix table of non-related individuals with PCA scores
    Remove related samples from PC relate from pruned MT and run PCA
    :param str pruned_mt_file: ld pruned MT file
    :param str pca_mt_file: PCA output MT file
    :param str mtdir: directory output matrix tables are written to
    :param str plotdir: directory output plots are written to
    :param str samples_to_remove_file: related samples to remove ht file
    '''
    print("Running population PCA")
    pruned_mt = hl.read_matrix_table(pruned_mt_file)
    print("Removing related samples")
    related_samples_to_remove = hl.read_table(samples_to_remove_file)
    pca_mt = pruned_mt.filter_cols(hl.is_defined(related_samples_to_remove[pruned_mt.col_key]), keep=False)
    variants, samples = pca_mt.count()
    print(f"{samples} samples after relatedness step.")

    plink_mt = pca_mt.annotate_cols(uid=pca_mt.s).key_cols_by('uid')
    plinkfile = mtdir + "mt_unrelated.plink"
    hl.export_plink(plink_mt, plinkfile,fam_id=plink_mt.uid, ind_id=plink_mt.uid)

    print("Running PCA")

    # TODO: DEBUG: k=20 didn't work for the small test dataset
    # Error summary: IllegalArgumentException: requirement failed: Requested k singular values but got k=20 and numCols=9.
    pca_evals, pca_scores, pca_loadings = hl.hwe_normalized_pca(pca_mt.GT, k=4, compute_loadings=True)
    pca_af_ht = pca_mt.annotate_rows(pca_af=hl.agg.mean(pca_mt.GT.n_alt_alleles()) / 2).rows()
    pca_loadings = pca_loadings.annotate(pca_af=pca_af_ht[pca_loadings.key].pca_af)
    pca_scores_file = mtdir + "mt_pca_scores.ht"#table of pca scores per sample
    pca_scores.write(pca_scores_file, overwrite=True)
    pca_loadings_file = mtdir + "mt_pca_loadings.ht"#table of SNPs, PCA loadings and AF
    pca_loadings.write(pca_loadings_file, overwrite=True)
    pca_mt = pca_mt.annotate_cols(scores=pca_scores[pca_mt.col_key].scores)
    pca_mt.write(pca_mt_file, overwrite=True)

    print("Plotting PC1 vs PC2")
    p = hl.plot.scatter(pca_mt.scores[0], pca_mt.scores[1], title='PCA', xlabel='PC1', ylabel='PC2')
    # TODO: change to os.path.join
    # TODO: add check if plotdir exists, create if not
    output_file(plotdir + "pca.html")
    save(p)



def main():
    #set up
    inputs = parse_config()
    mtdir = inputs['matrixtables_lustre_dir']
    plotdir = inputs['plots_dir_local']

    #initialise hail
    # TODO: clarify why these won't work for me
    # tmp_dir = "hdfs://spark-master:9820/"
    # tmp_dir = "hdfs://spark-master:9820/user/ubuntu/hail-tmp"

    tmp_dir = "file:///lustre/scratch126/dh24_test/tmp"

    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    #load input mt
    mtfile = mtdir + "mt_sex_annotated.mt"
    mt = hl.read_matrix_table(mtfile)

    #ld prune to get a table of variants which are not correlated
    pruned_mt_file = mtdir + "mt_ldpruned.mt"
    prune_mt(mt, pruned_mt_file)

    #run pcrelate
    relatedness_ht_file = mtdir + "mt_relatedness.ht"
    samples_to_remove_file = mtdir + "mt_related_samples_to_remove.ht"
    scores_file = mtdir + "mt_pruned.pca_scores.ht"
    run_pc_relate(pruned_mt_file, relatedness_ht_file, samples_to_remove_file, scores_file)

    #run PCA
    pca_mt_file = mtdir + "mt_pca.mt"
    run_population_pca(pruned_mt_file, pca_mt_file, mtdir, plotdir, samples_to_remove_file)


if __name__ == '__main__':
    main() 
