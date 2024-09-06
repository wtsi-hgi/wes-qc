# identify and prune related samples prior to PCA
# use mt with hard filters and sex annotation from 2-sample_qc/1-hard_filters_sex_annotation.py
import hail as hl
import os
import pyspark
from utils.utils import parse_config, subdict, multigetp, path_local, path_spark
import bokeh.plotting as bkplt


def prune_mt(mtin: hl.MatrixTable, config: dict) -> hl.MatrixTable:
    '''
    Splits multiallelic sites and runs ld pruning
    Filter to autosomes before LD pruning to decrease sample size - autosomes only wanted for later steps
    :param MatrixTable mtin: input MT containing variants to be pruned
    :param dict config:
    :return: Pruned MatrixTable
    :rtype: hl.MatrixTable

    `hl.ld_prune` returns a maximal subset of variants that are nearly uncorrelated within each window.
    Requires the dataset to contain only diploid genotype calls and no multiallelic variants.
    See also: https://hail.is/docs/0.2/methods/genetics.html#hail.methods.ld_prune

    ### Config fields
    step2.prune.pruned_mt_file : output path  
    step2.prune.r2 : float  

    ### Optional ld_prune arguments  
    step2.prune.ld_prune_args.bp_window_size : int optional  
    step2.prune.ld_prune_args.memory_per_core : int optional  
    step2.prune.ld_prune_args.keep_higher_maf : bool optional  
    step2.prune.ld_prune_args.block_size : int optional  
    '''
    conf = config['step2']['prune']
    ld_prune_args = conf['ld_prune_args']
    mtoutfile = conf['pruned_mt_file']

    print("Filtering to autosomes")
    mtin = mtin.filter_rows(mtin.locus.in_autosome())
    print("Splitting multiallelic sites")
    # mtin = hl.split_multi(mtin)
    mtin = hl.split_multi_hts(mtin)#this shouldn't do anything as only biallelic sites are used
    print("Performing LD pruning")
    pruned_ht = hl.ld_prune(mtin.GT, **ld_prune_args)
    pruned_mt = mtin.filter_rows(hl.is_defined(pruned_ht[mtin.row_key]))
    pruned_mt = pruned_mt.select_entries(GT=hl.unphased_diploid_gt_index_call(pruned_mt.GT.n_alt_alleles()))
    pruned_mt.write(path_spark(mtoutfile), overwrite=True) # output

    return pruned_mt


def run_pc_relate(pruned_mt: hl.MatrixTable, config: dict) -> hl.MatrixTable:
    '''
    Runs PC relate on pruned MT
    :param hl.MatrixTable pruned_mt: ld pruned MT file
    :param dict config:
    :return: mt with a column of too related samples to remove
    :rtype: hl.MatrixTable

    See https://hail.is/docs/0.2/methods/relatedness.html#hail.methods.pc_relate

    ### Config fields
    step2.pc_relate.relatedness_ht_file : output path : relatedness ht file  
    step2.pc_relate.samples_to_remove_file : output path : ht file with a list of samples to remove  
    step2.pc_relate.scores_file : output path: : ht file with hwe_normalized_pca scores  
    step2.pc_relate.pca_components : int : number of components for the hwe_normalized_pca   
    step2.pc_relate.relatedness_column : str : pc_relate output column to use when applying the relatedness_threshold  
    step2.pc_relate.relatedness_threshold : float : samples with relatedness_column greater than relatedness_threshold are to be removed  

    ### Optional config fields
    step2.pc_relate.pc_relate_args.min_individual_maf : float  
    step2.pc_relate.pc_relate_args.block_size : float  
    step2.pc_relate.pc_relate_args.min_kinship : float  
    step2.pc_relate.pc_relate_args.statistics : str  
    step2.pc_relate.pc_relate_args.k : int  
    step2.pc_relate.pc_relate_args.include_self_kinship : bool  
    '''
    conf = config['step2']['prune_pc_relate']
    pc_relate_kwargs = conf['pc_relate_args']

    print("Running PC relate")

    # TODO: DEBUG: k=10 didn't work for the small test dataset, throwing an error:
    # "Error summary: IllegalArgumentException: requirement failed: Dimension mismatch!: _a.cols == _b.rows (11 != 10)"
    eig, scores, _ = hl.hwe_normalized_pca(pruned_mt.GT, k=conf['pca_components'], compute_loadings=False)
    scores.write(path_spark(conf['scores_file']), overwrite=True) # output

    print("Calculating relatedness")
    relatedness_ht = hl.pc_relate(pruned_mt.GT, scores_expr=scores[pruned_mt.col_key].scores, **pc_relate_kwargs)
    relatedness_ht.write(path_spark(conf['relatedness_ht_file']), overwrite=True) # output
    
    # prune individuals to be left with unrelated - creates a table containing one column - samples to remove
    rel_col = conf['relatedness_column']
    rel_threshold = conf['relatedness_threshold']

    pairs = relatedness_ht.filter(relatedness_ht[rel_col] > rel_threshold)
    related_samples_to_remove = hl.maximal_independent_set(pairs.i, pairs.j, keep=False)
    related_samples_to_remove.write(path_spark(conf['samples_to_remove_file']), overwrite=True) # output

    return related_samples_to_remove


# How is this step different from 2-sample_qc/3-population_pca_prediction.py / run_pca ?
# Only PCA plotting?
def run_population_pca(pruned_mt: hl.MatrixTable, samples_to_remove: hl.Table, config: dict) -> hl.MatrixTable:
    '''
    Runs PCA and creates a matrix table of non-related individuals with PCA scores
    Remove related samples from PC relate from pruned MT and run PCA
    :param hl.MatrixTable pruned_mt: ld pruned MT
    :param hl.Table samples_to_remove_file: a HT with a single column of related samples to remove 
    :param dict config:
    :param str pca_mt_file: PCA output MT file path TODO move to config
    :return: TODO
    :rtype: hl.MatrixTable

    ### Config fields
    step2.population_pca.plink_outfile : output path : TODO  
    step2.population_pca.pca_components : int : TODO  
    step2.population_pca.pca_scores_file : output path : TODO  
    step2.population_pca.pca_loadings_file : output path : TODO  
    step2.population_pca.pca_mt_file : output path : TODO  
    step2.population_pca.plot_outfile : output path : Bokeh PCA scatterplot  

    '''
    conf = config['step2']['prune_plot_pca']
    plinkfile = conf['plink_outfile']
    pca_scores_file = conf['pca_scores_file']
    pca_loadings_file = conf['pca_loadings_file']
    pca_mt_file = conf['pca_mt_file']

    print("Running population PCA")
    print("Removing related samples")
    pca_mt = pruned_mt.filter_cols(hl.is_defined(samples_to_remove[pruned_mt.col_key]), keep=False)
    variants, samples = pca_mt.count()
    print(f"{samples} samples after relatedness step.")

    plink_mt = pca_mt.annotate_cols(uid=pca_mt.s).key_cols_by('uid')
    hl.export_plink(dataset=plink_mt, output=path_spark(plinkfile), fam_id=plink_mt.uid, ind_id=plink_mt.uid) # output

    print("Running PCA")

    # TODO: DEBUG: k=20 didn't work for the small test dataset
    # Error summary: IllegalArgumentException: requirement failed: Requested k singular values but got k=20 and numCols=9.
    pca_evals, pca_scores, pca_loadings = hl.hwe_normalized_pca(pca_mt.GT, k=conf['pca_components'], compute_loadings=True)
    pca_af_ht = pca_mt.annotate_rows(pca_af=hl.agg.mean(pca_mt.GT.n_alt_alleles()) / 2).rows()
    pca_loadings = pca_loadings.annotate(pca_af=pca_af_ht[pca_loadings.key].pca_af)
    pca_scores.write(path_spark(pca_scores_file), overwrite=True) # output
    pca_loadings.write(path_spark(pca_loadings_file), overwrite=True) # output
    pca_mt = pca_mt.annotate_cols(scores=pca_scores[pca_mt.col_key].scores)
    pca_mt.write(path_spark(pca_mt_file), overwrite=True) # output

    print("Plotting PC1 vs PC2")
    plot_outfile = path_local(conf['plot_outfile'])
    p = hl.plot.scatter(pca_mt.scores[0], pca_mt.scores[1], title='PCA', xlabel='PC1', ylabel='PC2')
    print(f'saving plot to {plot_outfile} DEBUG') # DEBUG
    bkplt.output_file(plot_outfile)
    bkplt.save(p)

    return pca_mt



def main():
    #set up
    config = parse_config()

    #initialise hailS
    tmp_dir = config['general']['tmp_dir']
    # sc = pyspark.SparkContext()
    sc = pyspark.SparkContext.getOrCreate()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38", idempotent=True)

    #ensure plotting directory exists
    pltdir = path_local(config['general']['plots_dir'])
    if not os.path.exists(pltdir):
        os.makedirs(pltdir)

    #load input mt
    # mt_infile = os.path.join(mtdir, "mt_sex_annotated.mt")
    mt_infile = config['step2']['impute_sex']['sex_mt_outfile']
    mt = hl.read_matrix_table(path_spark(mt_infile))

    #ld prune to get a table of variants which are not correlated
    # pruned_mt_file = os.path.join(mtdir, "mt_ldpruned.mt")
    pruned_mt = prune_mt(mt, config)

    #run pcrelate
    # relatedness_ht_file = os.path.join(mtdir, "mt_relatedness.ht")
    # samples_to_remove_file = os.path.join(mtdir, "mt_related_samples_to_remove.ht")
    # scores_file = os.path.join(mtdir, "mt_pruned.pca_scores.ht")
    related_samples_to_remove_ht = run_pc_relate(pruned_mt, config)

    #run PCA
    population_pca_mt = run_population_pca(pruned_mt, related_samples_to_remove_ht, config)


if __name__ == '__main__':
    main() 
