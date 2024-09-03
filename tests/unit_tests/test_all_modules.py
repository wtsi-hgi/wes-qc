import os
import re
import gzip
import unittest
import importlib
import hail as hl
import hailtop.fs as hfs
import shutil as sh
from pyspark import SparkContext
from utils.config import parse_config, path_local, path_spark, getp, _expand_cvars_recursively
import subprocess


def compare_structs(struct1, struct2):
    """
    Compare each field in two struct entries.
    
    Returns:
    bool: True if entries are equal, False otherwise.
    """
    fields = list(struct1.dtype.fields)

    # if not fields other than key
    if len(fields) == 0:
        return hl.bool(True)
    
    # Check for NA fields
    def check_field_equal(field):
        # expr that field is NA in both entries (MT) or rows (HT)
        field_not_defined = (~hl.is_defined(struct1[field]) & ~hl.is_defined(struct2[field]))
        # expr that field is not NA in both entries (MT) or rows (HT)
        field_defined = hl.is_defined(struct1[field]) & hl.is_defined(struct2[field])
        # different checks for differenet field dtypes
        if struct1[field].dtype == hl.dtype('float'):
            expr = (field_defined & hl.approx_equal(struct1[field], struct2[field])) | field_not_defined
        elif struct1[field].dtype == hl.dtype('array<float>'):
            # element_wise_expr = hl.map(lambda a, b: (hl.is_defined(a) & hl.is_defined(b) & hl.approx_equal(a, b)) | (~hl.is_defined(a) & ~hl.is_defined(b)), struct1[field], struct2[field])

            element_wise_expr = hl.map(lambda a, b: (hl.is_defined(a) &
                                                     hl.is_defined(b) &
                                                     hl.approx_equal(a, b)) | 
                                                     (~hl.is_defined(a) & 
                                                      ~hl.is_defined(b)), 
                                       struct1[field], struct2[field])
            expr = hl.all(element_wise_expr)
        else:
            # TODO: test for other dtypes
            expr = (field_defined & (struct1[field] == struct2[field])) | field_not_defined
        
        return expr

    comparisons = [check_field_equal(field) for field in fields]
    return hl.all(lambda x: x, comparisons)

def compare_entries_to_other_mt(this_entry, row_key, col_key, other_mt):
    """
    Compare entries of two MatrixTables.
    
    Returns:
    bool: True if entries are equal, False otherwise.
    """
    other_entry = other_mt[row_key, col_key]
    return ~compare_structs(this_entry, other_entry)

def compare_plinks(bed1_path: str, bim1_path: str, fam1_path: str, 
                   bed2_path: str, bim2_path: str, fam2_path: str) -> bool:
    """
    Compare the contents of two Hail plink files.

    Returns:
    bool: True if the plinks are equal, False otherwise.
    """
    # Load the plinks
    mt1 = hl.import_plink(bed=bed1_path, bim=bim1_path, fam=fam1_path, reference_genome='GRCh38')
    mt2 = hl.import_plink(bed=bed2_path, bim=bim2_path, fam=fam2_path, reference_genome='GRCh38')

    # Ensure the schemas are the same
    if mt1.row.dtype != mt2.row.dtype or mt1.col.dtype != mt2.col.dtype or mt1.entry.dtype != mt2.entry.dtype:
        print('MatrixTable schemas do not match')
        return False

    # Align row and column keys
    mt2 = mt2.key_rows_by(*mt1.row_key)  # Reorder row keys to match mt1
    mt2 = mt2.key_cols_by(*mt1.col_key)  # Reorder column keys to match mt1
    
    filtered_mt = mt1.filter_entries(compare_entries_to_other_mt(mt1.entry, mt1.row_key, mt1.col_key, mt2))

    # Check if there are any differing entries
    num_differences = filtered_mt.entries().count()
    if num_differences == 0:
        print('MatrixTables are equal') # DEBUG
        return True
    else:
        print(f'MatrixTables are not equal: {num_differences} differing entries found')
        return False

def compare_matrixtables(mt1_path: str, mt2_path: str) -> bool:
    """
    Compare the contents of two Hail MatrixTables.

    Parameters:
    mt1_path (str): Path to the first MatrixTable.
    mt2_path (str): Path to the second MatrixTable.

    Returns:
    bool: True if the MatrixTables are equal, False otherwise.
    """
    # Load the MatrixTables
    mt1 = hl.read_matrix_table(mt1_path)
    mt2 = hl.read_matrix_table(mt2_path)

    # Ensure the schemas are the same
    if mt1.row.dtype != mt2.row.dtype or mt1.col.dtype != mt2.col.dtype or mt1.entry.dtype != mt2.entry.dtype:
        print('MatrixTable schemas do not match')
        return False

    # Align row and column keys
    mt2 = mt2.key_rows_by(*mt1.row_key)  # Reorder row keys to match mt1
    mt2 = mt2.key_cols_by(*mt1.col_key)  # Reorder column keys to match mt1
    
    filtered_mt = mt1.filter_entries(compare_entries_to_other_mt(mt1.entry, mt1.row_key, mt1.col_key, mt2))

    # Check if there are any differing entries
    num_differences = filtered_mt.entries().count()
    if num_differences == 0:
        print('MatrixTables are equal') # DEBUG
        return True
    else:
        print(f'MatrixTables are not equal: {num_differences} differing entries found')
        return False

def compare_entries_to_other_ht(this_row_value, row_key, other_ht):
    """
    Compare rows of two Tables.
    
    Returns:
    bool: True if entries are equal, False otherwise.
    """
    other_row_value = other_ht[row_key]
    return ~compare_structs(this_row_value, other_row_value)

def compare_tables(ht1_path: str, ht2_path: str) -> bool:
    """
    Compare the contents of two Hail Tables.

    Parameters:
    ht1_path (str): Path to the first Table.
    ht2_path (str): Path to the second Table.

    Returns:
    bool: True if the Tables are equal, False otherwise.
    """
    # Load the Tables
    ht1 = hl.read_table(ht1_path)
    ht2 = hl.read_table(ht2_path)

    # Ensure the schemas are the same
    if ht1.row.dtype != ht2.row.dtype:
        print("Table schemas do not match")
        return False

    # Align row keys
    ht2 = ht2.key_by(*ht1.key)  # Reorder row keys to match ht1

    rows_differ = ht1.filter(compare_entries_to_other_ht(ht1.row_value, ht1.key, ht2))

    # Check if there are any differing rows
    num_differences = rows_differ.count()

    if num_differences == 0:
        print(f'Tables {ht1_path}\n{ht2_path}\nare equal') # DEBUG
        return True
    else:
        print(f'Tables {ht1_path}\n{ht2_path}\nare not equal: {num_differences} differing rows found') # DEBUG
        return False

def compare_txts(path_1: str, path_2: str, replace_strings: list[list[str, str]] = None) -> bool:
    """
    Compare contents of two files based on their raw string content.

    Returns:
    bool: True if files' contents are equal, False otherwise.
    """
    with open(path_1, 'r') as f_1:
        contents_1 = f_1.read()

    with open(path_2, 'r') as f_2:
        contents_2 = f_2.read()

    if replace_strings:
        for pattern, substitute in replace_strings:
            contents_1 = re.sub(pattern, substitute, contents_1)
            contents_2 = re.sub(pattern, substitute, contents_2)

    if contents_1 == contents_2:
        print(f'Files {path_1}\n{path_2}\n are equal') # DEBUG
    else:
        print(f'Files {path_1}\n{path_2}\n are not equal')

    return contents_1 == contents_2

def compare_bgzed_txts(path_1: str, path_2: str, replace_strings: list[list[str, str]] = None) -> bool:
    """
    Compare contents of two bgzed txt files. 
    
    Parameters:
    replace_strings: list[list[str, str]] - regex and str to replace

    Returns:
    bool: True if files' contents are equal, False otherwise.
    """
    with gzip.open(path_1, 'r') as f_1:
        contents_1 = f_1.read()

    with gzip.open(path_2, 'r') as f_2:
        contents_2 = f_2.read()
    
    if replace_strings:
        for pattern, substitute in replace_strings:
            contents_1 = re.sub(pattern, substitute, contents_1)
            contents_2 = re.sub(pattern, substitute, contents_2)

    if contents_1 == contents_2:
        print(f'Files {path_1}\n{path_2}\n are equal') # DEBUG
    else:
        print(f'Files {path_1}\n{path_2}\n are not equal')

    return contents_1 == contents_2

# ensure sequential test execution order
# unittest.TestLoader.sortTestMethodsUsing = None

# DEBUG: might be useful https://spark.apache.org/docs/latest/api/python/getting_started/testing_pyspark.html
# TODO: clean up hail logs

# /path/to/wes_qc must be in PYTHONPATH
qc_step_1_1 = importlib.import_module("1-import_data.1-import_gatk_vcfs_to_hail")
qc_step_2_1 = importlib.import_module("2-sample_qc.1-hard_filters_sex_annotation")
qc_step_2_2 = importlib.import_module("2-sample_qc.2-prune_related_samples")
qc_step_2_3 = importlib.import_module("2-sample_qc.3-population_pca_prediction")
qc_step_2_4 = importlib.import_module("2-sample_qc.4-find_population_outliers")
qc_step_2_5 = importlib.import_module("2-sample_qc.5-filter_fail_sample_qc")


class HailTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # parent outdir
        cls.test_outdir_path = os.path.dirname(os.path.realpath(__file__))

        # tmp dir for Hail
        cls.tmp_dir = f"file://{os.path.join(cls.test_outdir_path, 'tmp_test')}"

        sc = SparkContext.getOrCreate()
        # hadoop_config = sc._jsc.hadoopConfiguration() # DEBUG not used, why?

        # initialise Hail in local mode to run the test even without cluster
        hl.init(sc=sc, app_name="HailTest", master="local[*]",
                default_reference="GRCh38", tmp_dir=cls.tmp_dir)

    @classmethod
    def tearDownClass(cls):
        # clean up tmp_dir
        hfs.rmtree(cls.tmp_dir)
        hl.stop()


class TestQCSteps(HailTestCase):
    @classmethod
    def setUpClass(cls):
        # set up Hail from the base class
        super(TestQCSteps, cls).setUpClass()

        # define parameters needed for functions to be tested

        # path to dir with the test data in the repo
        cls.test_dataset_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        cls.ref_dataset_path = f"file://{os.path.join(os.path.dirname(os.path.realpath(__file__)), 'reference_output_data')}"
        refdir = cls.ref_dataset_path.removeprefix('file://')
        os.makedirs(refdir, exist_ok=True)
        subprocess.run(['s3cmd', 'get', '-r', '--skip-existing', 's3://wes-qc-data/unit_tests/reference_output_data/', refdir]) #BEEP

        cls.test_resourcedir = f"file://{os.path.join(cls.test_dataset_path, 'resources')}"
        resdir = cls.test_resourcedir.removeprefix('file://')
        os.makedirs(resdir, exist_ok=True)
        subprocess.run(['s3cmd', 'get', '-r', '--skip-existing', 's3://wes-qc-data/resources/', resdir]) #BEEP
        cls.onekg_resdir = f"file://{os.path.join(cls.test_dataset_path, 'resources', 'mini_1000G')}"


        # TODO: switch to config rendering as in integration tests?

        # TODO: move below (near the config setup) after refactoring other steps
        # ===== QC General Params ===== #
        cls.mtdir = f"file://{os.path.join(cls.test_outdir_path, 'matrixtables_test')}"
        cls.annotdir = f"file://{os.path.join(cls.test_outdir_path, 'annotations_test')}"
        cls.plots_dir = f"file://{os.path.join(cls.test_outdir_path, 'plots_test')}"
        os.makedirs(cls.plots_dir.removeprefix('file://'), exist_ok=True)

        # mk43
        # TODO: use a template engine to create a yaml
        # TODO: cleanup this function
        # TODO: separate this functions into steps
        
        # add general params into the config object
        config = dict()
        config['general'] = dict(
            tmp_dir = cls.tmp_dir,
            annotation_dir = cls.annotdir,
            matrixtables_dir = cls.mtdir,
            resource_dir = cls.test_resourcedir,
            plots_dir = cls.plots_dir,
            onekg_resource_dir = cls.onekg_resdir
        )

        # ===== QC Step 1.1 ===== #
        # inputs
        cls.import_vcf_dir = f"file://{os.path.join(cls.test_dataset_path, 'control_set_small')}"
        vcfdir = cls.import_vcf_dir.removeprefix('file://')
        os.makedirs(vcfdir, exist_ok=True)
        subprocess.run(['s3cmd', 'get', '-r', '--skip-existing', 's3://wes-qc-data/control_set_small/', vcfdir]) #BEEP
        
        cls.vcf_header = '' # not available for test dataset. TODO: create one to ensure test completeness
        # outputs
        # TODO: make output filename variable and refactor the test
        cls.output_mt_path = os.path.join(cls.mtdir, 'gatk_unprocessed.mt') # DEBUG: keep for testing new config format
        # reference outputs
        cls.ref_mt_file = os.path.join(cls.ref_dataset_path, 'matrixtables', 'gatk_unprocessed.mt')

        # add step 1 parameters into the config obj
        config['step1'] = dict(
            gatk_vcf_header_infile = cls.vcf_header,
            gatk_vcf_indir = cls.import_vcf_dir,
            gatk_mt_outfile = '{mtdir}/gatk_unprocessed.mt',
        )

        # ===== QC Step 2.1 ===== #
        # apply_hard_filters_outputs()
        # reference inputs
        cls.ref_unfiltered_mt_path = os.path.join(cls.ref_dataset_path, 'matrixtables', 'gatk_unprocessed.mt')
        # outputs
        # TODO: make output filename variable and refactor the test
        cls.output_filtered_mt_path = os.path.join(cls.mtdir, 'mt_hard_filters_annotated.mt') # DEBUG: keep for testing new config format
        # reference outputs
        cls.ref_output_filtered_mt_path = os.path.join(cls.ref_dataset_path, 'matrixtables', 'mt_hard_filters_annotated.mt')

        # impute_sex()
        # reference inputs: `cls.ref_output_filtered_mt_path`
        # outputs
        cls.output_sex_annotated_path = os.path.join(cls.annotdir, 'sex_annotated.sex_check.txt.bgz') # DEBUG: keep for testing new config format
        cls.output_sex_mt_path = os.path.join(cls.mtdir, 'mt_sex_annotated.mt') # DEBUG: keep for testing new config format
        # reference outputs
        cls.ref_output_sex_annotated_path = os.path.join(cls.ref_dataset_path, 'annotations', 'sex_annotated.sex_check.txt.bgz')
        cls.ref_output_sex_mt_path = os.path.join(cls.ref_dataset_path, 'matrixtables', 'mt_sex_annotated.mt')

        # identify_inconsistencies()
        # reference inputs: `cls.ref_output_sex_mt_path`
        # outputs
        cls.conflicting_sex_path = os.path.join(cls.annotdir, 'conflicting_sex.txt.bgz') # DEBUG: keep for testing new config format
        cls.sex_annotation_f_stat_outliers_path = os.path.join(cls.annotdir, 'sex_annotation_f_stat_outliers.txt.bgz') # DEBUG: keep for testing new config format
        # reference outputs
        cls.ref_conflicting_sex_path = os.path.join(cls.ref_dataset_path, 'annotations', 'conflicting_sex.txt.bgz')
        cls.ref_sex_annotation_f_stat_outliers_path = os.path.join(cls.ref_dataset_path, 'annotations', 'sex_annotation_f_stat_outliers.txt.bgz')

        # add step 2.1 parameters into the config obj
        # TODO: create variables for threshold values
        config['step2'] = dict()
        config['step2']['sex_annotation_hard_filters'] = dict(
            filtered_mt_outfile = '{mtdir}/mt_hard_filters_annotated.mt',
            n_alt_alleles_threshold = 0.001,
            defined_gt_frac_threshold = 0.99
        )
        config['step2']['impute_sex'] = {
            'sex_ht_outfile': '{anndir}/sex_annotated.sex_check.txt.bgz',
            'sex_mt_outfile': '{mtdir}/mt_sex_annotated.mt',
            'female_threshold': 0.5,
            'male_threshold': 0.8,
            'aaf_threshold': 0.05
        }
        config['step2']['sex_inconsistencies'] = {
            'sex_metadata_file': '{resdir}/mlwh_sample_and_sex.txt',
            'conflicting_sex_report_file': '{anndir}/conflicting_sex.txt.bgz',
            'fstat_outliers_report_file': '{anndir}/sex_annotation_f_stat_outliers.txt.bgz',
            'fstat_low': 0.2,
            'fstat_high': 0.8
        }

        # ===== QC Step 2.2 ===== #
        # prune_mt()
        # reference inputs: `cls.ref_output_sex_mt_path`
        # outputs
        cls.mt_ldpruned_path = os.path.join(cls.mtdir, 'mt_ldpruned.mt')
        # reference outputs
        cls.ref_mt_ldpruned_path = os.path.join(cls.ref_dataset_path, 'matrixtables', 'mt_ldpruned.mt')

        # run_pc_relate()
        # reference inputs: `cls.ref_mt_ldpruned_path`
        # outputs
        cls.relatedness_ht_path = os.path.join(cls.mtdir, 'mt_relatedness.ht')
        cls._path_path = os.path.join(cls.mtdir, 'mt_related_samples_to_remove.ht')
        cls.scores_path = os.path.join(cls.mtdir, 'mt_pruned.pca_scores.ht')
        # reference outputs
        cls.ref_relatedness_ht_path = os.path.join(cls.ref_dataset_path, 'matrixtables', 'mt_relatedness.ht')
        cls.ref_samples_to_remove_path = os.path.join(cls.ref_dataset_path, 'matrixtables', 'mt_related_samples_to_remove.ht')
        cls.ref_scores_path = os.path.join(cls.ref_dataset_path, 'matrixtables', 'mt_pruned.pca_scores.ht')

        # run_population_pca()
        # reference inputs: `cls.ref_mt_ldpruned_path`, `cls.ref_samples_to_remove_path`
        # outputs
        cls.plotdir = os.path.join(cls.test_outdir_path, 'plots_test')
        cls.plink_path = os.path.join(cls.mtdir, 'mt_unrelated.plink')
        cls.mt_pca_scores_path = os.path.join(cls.mtdir, 'mt_pca_scores.ht')
        cls.mt_pca_loadings_path = os.path.join(cls.mtdir, 'mt_pca_loadings.ht')
        cls.pca_output_path = os.path.join(cls.plotdir, 'pca.html')
        cls.pca_mt_path = os.path.join(cls.mtdir, 'mt_pca.mt')
        # reference outputs
        cls.ref_plink_path = os.path.join(cls.ref_dataset_path, 'matrixtables', 'mt_unrelated.plink')
        cls.ref_mt_pca_scores_path = os.path.join(cls.ref_dataset_path, 'matrixtables', 'mt_pca_scores.ht')
        cls.ref_mt_pca_loadings_path = os.path.join(cls.ref_dataset_path, 'matrixtables', 'mt_pca_loadings.ht')
        cls.ref_pca_output_path = path_local(os.path.join(cls.ref_dataset_path, 'plots', 'pca.html'))
        cls.ref_pca_mt_path = os.path.join(cls.ref_dataset_path, 'matrixtables', 'mt_pca.mt')

        config['step2']['prune'] = {
            'pruned_mt_outfile': '{mtdir}/mt_ldpruned.mt',
            'ld_prune_args': {'r2': 0.2}
        }

        config['step2']['prune_pc_relate'] = {
            'relatedness_ht_file' : '{mtdir}/mt_relatedness.ht',
            'samples_to_remove_file' : '{mtdir}/mt_related_samples_to_remove.ht',
            'scores_file': '{mtdir}/mt_pruned.pca_scores.ht',
            'pca_components': 3,
            'pc_relate_args': {'min_individual_maf': 0.05,
                               'block_size' : 4096,
                               'min_kinship': 0.05,
                               'statistics': 'kin2',
                               # k:,
                               # include_self_kinship:
                               },         
            'relatedness_column' : 'kin',
            'relatedness_threshold': 0.125
        }

        config['step2']['prune_plot_pca'] = {
            'plink_outfile' : '{mtdir}/mt_unrelated.plink',
            'pca_components' : 4,
            'pca_scores_file' : '{mtdir}/mt_pca_scores.ht',
            'pca_loadings_file' : '{mtdir}/mt_pca_loadings.ht',
            'pca_mt_file' : '{mtdir}/mt_pca.mt',
            'plot_outfile' : '{pltdir}/pca.html'
        }

        # ===== QC Step 2.3 ===== #
        # create_1kg_mt()
        # reference inputs: no
        # outputs
        cls.kg_wes_regions = os.path.join(cls.mtdir, 'kg_wes_regions.mt')
        # reference outputs
        cls.ref_kg_wes_regions = os.path.join(cls.ref_dataset_path, 'matrixtables', 'kg_wes_regions.mt')

        # merge_with_1kg()
        # reference inputs: `cls.ref_mt_ldpruned_path`
        cls.ref_mtdir = os.path.join(cls.ref_dataset_path, 'matrixtables')
        # outputs
        cls.merged_mt_path = os.path.join(cls.mtdir, 'merged_with_1kg.mt')
        # reference outputs
        cls.ref_merged_mt_path = os.path.join(cls.ref_mtdir, "merged_with_1kg.mt")

        # annotate_and_filter()
        # reference inputs: `cls.ref_merged_mt_path`, `cls.test_resourcedir`
        # outputs
        cls.filtered_mt_path = os.path.join(cls.mtdir, 'merged_with_1kg_filtered.mt')
        # reference outputs
        cls.ref_filtered_mt_path = os.path.join(cls.ref_mtdir, "merged_with_1kg_filtered.mt")

        # run_pca()
        # reference inputs: `cls.ref_filtered_mt_path`
        # outputs
        cls.pca_scores_path = os.path.join(cls.mtdir, 'pca_scores_after_pruning.ht')
        cls.pca_loadings_path = os.path.join(cls.mtdir, 'pca_loadings_after_pruning.ht')
        cls.pca_evals_path = path_local(os.path.join(cls.mtdir, 'pca_evals_after_pruning.txt'))
        # reference outputs
        cls.ref_pca_scores_path = os.path.join(cls.ref_mtdir, 'pca_scores_after_pruning.ht')
        cls.ref_pca_loadings_path = os.path.join(cls.ref_mtdir, 'pca_loadings_after_pruning.ht')
        cls.ref_pca_evals_path = path_local(os.path.join(cls.ref_mtdir, 'pca_evals_after_pruning.txt'))

        # predict_pops()
        # reference inputs: `cls.ref_pca_scores_path`
        # outputs
        cls.pop_ht_file = os.path.join(cls.mtdir, 'pop_assignments.ht')
        cls.pop_ht_tsv = os.path.join(cls.mtdir, 'pop_assignments.tsv')
        # reference outputs
        cls.ref_pop_ht_file = os.path.join(cls.ref_mtdir, 'pop_assignments.ht')
        cls.ref_pop_ht_tsv = os.path.join(cls.ref_mtdir, 'pop_assignments.tsv')

        config['step2']['create_1kg_mt'] = {
            'indir': '{1kg_resdir}',
            'vcfheader': '{1kg_resdir}/header_20201028.txt',
            'mt_out_file' : '{mtdir}/kg_wes_regions.mt'
        }

        config['step2']['merge_with_1_kg'] = {
            'kg_mt_file' : '{mtdir}/kg_wes_regions.mt',
            'merged_mt_outfile' : '{mtdir}/merged_with_1kg.mt'
        }

        config['step2']['annotate_and_filter'] = {
            'pops_file' : '{resdir}/igsr_samples.tsv',
            'call_rate' : 0.99,
            'AF' : 0.05,
            'p_value_hwe' : 0.00005,
            'long_range_ld_file' : '{resdir}/long_ld_regions.hg38.bed',
            'filtered_mt_outfile' : '{mtdir}/merged_with_1kg_filtered.mt'
        }

        config['step2']['run_pca'] = {
            'pca_scores_outfile' : '{mtdir}/pca_scores_after_pruning.ht',
            'pca_loadings_outfile' : '{mtdir}/pca_loadings_after_pruning.ht',
            'pca_evals_outfile' : '{mtdir}/pca_evals_after_pruning.txt'
        }
        
        config['step2']['predict_pops'] = {
            'pca_scores_file' : '{mtdir}/pca_scores_after_pruning.ht',
            'gnomad_pc_n_estimators' : 100,
            'gnomad_prop_train' : 0.8,
            'gnomad_min_prob' : 0.5,
            'pop_ht_outfile' : '{mtdir}/pop_assignments.ht',
            'pop_ht_outtsv' : '{mtdir}/pop_assignments.tsv'}

        # ===== QC Step 2.4 ===== #
        # annotate_mt()
        # reference inputs: `cls.ref_pop_ht_file`, `cls.ref_pop_ht_tsv`, `cls.ref_mt_file`
        # outputs
        cls.annotated_mt_file = os.path.join(cls.mtdir, 'gatk_unprocessed_with_pop.mt')
        # reference outputs
        cls.ref_annotated_mt_file = os.path.join(cls.ref_mtdir, 'gatk_unprocessed_with_pop.mt')

        # stratified_sample_qc()
        # reference inputs: `cls.ref_annotated_mt_file`, `cls.ref_annotdir`
        cls.ref_annotdir = os.path.join(cls.ref_dataset_path, 'annotations')
        # outputs
        cls.mt_qc_outfile = os.path.join(cls.mtdir, 'mt_pops_sampleqc.mt')
        cls.ht_qc_cols_outfile = os.path.join(cls.mtdir, 'mt_pops_sampleqc.ht')
        cls.qc_filter_file = os.path.join(cls.mtdir, 'mt_pops_QC_filters.ht')
        cls.output_text_file = os.path.join(cls.annotdir, 'sample_qc_by_pop.tsv.bgz')
        cls.output_globals_json = os.path.join(cls.annotdir, 'sample_qc_by_pop.globals.json')

        # reference outputs
        cls.ref_mt_qc_outfile = os.path.join(cls.ref_mtdir, 'mt_pops_sampleqc.mt')
        cls.ref_ht_qc_cols_outfile = os.path.join(cls.ref_mtdir, 'mt_pops_sampleqc.ht')
        cls.ref_qc_filter_file = os.path.join(cls.ref_mtdir, 'mt_pops_QC_filters.ht')
        cls.ref_output_text_file = os.path.join(cls.ref_annotdir, 'sample_qc_by_pop.tsv.bgz')
        cls.ref_output_globals_json = os.path.join(cls.ref_annotdir, 'sample_qc_by_pop.globals.json')
        
        config['step2']['annotate_with_pop'] = {
            'annotated_mt_file' : '{mtdir}/gatk_unprocessed_with_pop.mt'
        }

        config['step2']['stratified_sample_qc'] = {
            'mt_qc_outfile': '{mtdir}/mt_pops_sampleqc.mt',
            'ht_qc_cols_outfile' : '{mtdir}/mt_pops_sampleqc.ht',
            'qc_filter_file' : '{mtdir}/mt_pops_QC_filters.ht',
            'min_depth' : 20,
            'min_genotype_quality' : 20,
            'min_vaf' : 0.25,
            'output_text_file' : '{anndir}/sample_qc_by_pop.tsv.bgz',
            'output_globals_json_file' : '{anndir}/sample_qc_by_pop.globals.json'
        }

        # ===== QC Step 2.5 ===== #
        # remove_sample_qc_fails()
        # reference inputs: `cls.ref_annotated_mt_file`, cls.ref_qc_filter_file`
        # outputs
        cls.sample_qc_filtered_mt_file = os.path.join(cls.mtdir, 'mt_pops_QC_filters_after_sample_qc.mt')
        cls.samples_failing_qc_file = os.path.join(cls.annotdir, 'samples_failing_qc.tsv.bgz')
        # reference outputs
        cls.ref_sample_qc_filtered_mt_file = os.path.join(cls.ref_mtdir, 'mt_pops_QC_filters_after_sample_qc.mt')
        cls.ref_samples_failing_qc_file = os.path.join(cls.ref_annotdir, 'samples_failing_qc.tsv.bgz')

        config['step2']['remove_sample_qc_fails'] = {
            'samples_failing_qc_file' : '{anndir}/samples_failing_qc.tsv.bgz',
            'sample_qc_filtered_mt_file' : '{mtdir}/mt_pops_QC_filters_after_sample_qc.mt'
        }

        cls.config = _expand_cvars_recursively(config, config)
        
    # QC Step 1.1
    def test_1_1_load_vcfs_to_mt(self):
        # run function to test
        qc_step_1_1.load_vcfs_to_mt(self.config)
        
        # compare the output to reference
        output_mt_path = self.config['step1']['gatk_mt_outfile']
        self.assertEqual(output_mt_path, self.output_mt_path) # DEBUG

        outputs_are_identical = compare_matrixtables(self.ref_mt_file, output_mt_path)
        self.assertTrue(outputs_are_identical)

    # QC Step 2.1 apply_hard_filters
    def test_2_1_1_apply_hard_filters(self):
        # read reference output of the step 1.1
        ref_mt_unfiltered = hl.read_matrix_table(self.ref_unfiltered_mt_path)

        # run function to test
        mt_filtered = qc_step_2_1.apply_hard_filters(ref_mt_unfiltered, self.config)

        
        # compare the output to reference
        output_filtered_mt_path = self.config['step2']['sex_annotation_hard_filters']['filtered_mt_outfile']
        self.assertEqual(output_filtered_mt_path, self.output_filtered_mt_path) # DEBUG
        
        outputs_are_identical = compare_matrixtables(self.ref_output_filtered_mt_path, output_filtered_mt_path)
        self.assertTrue(outputs_are_identical)

    # QC Step 2.1 impute sex
    def test_2_1_2_impute_sex(self):
        # read reference output of the Step 2.1 apply_hard_filters()
        ref_mt_filtered = hl.read_matrix_table(self.ref_output_filtered_mt_path)

        # run function to test
        sex_mt = qc_step_2_1.impute_sex(ref_mt_filtered, self.config)

        # compare the outputs to reference
        output_sex_annotated_path = self.config['step2']['impute_sex']['sex_ht_outfile']
        output_sex_mt_path = self.config['step2']['impute_sex']['sex_mt_outfile']
        self.assertEqual(output_sex_annotated_path, self.output_sex_annotated_path) # DEBUG
        self.assertEqual(output_sex_mt_path, self.output_sex_mt_path) # DEBUG

        self.assertEqual(path_local(output_sex_annotated_path), path_local(output_sex_annotated_path))

        output_mts_are_identical = compare_matrixtables(self.ref_output_sex_mt_path, output_sex_mt_path)
        output_txts_are_identical = compare_bgzed_txts(path_local(self.ref_output_sex_annotated_path),
                                                       path_local(output_sex_annotated_path))

        self.assertTrue(output_mts_are_identical and output_txts_are_identical)

    def test_2_1_3_identify_inconsistencies(self):
        # read reference output of the Step 2.1 impute_sex()
        ref_output_sex_mt = hl.read_matrix_table(self.ref_output_sex_mt_path)

        # run function to test
        qc_step_2_1.identify_inconsistencies(ref_output_sex_mt, self.config)

        # compare outputs to reference
        conflicting_sex_path = self.config['step2']['sex_inconsistencies']['conflicting_sex_report_file']
        sex_annotation_f_stat_outliers_path = self.config['step2']['sex_inconsistencies']['fstat_outliers_report_file']

        self.assertEqual(conflicting_sex_path, self.conflicting_sex_path) # DEBUG
        self.assertEqual(sex_annotation_f_stat_outliers_path, self.sex_annotation_f_stat_outliers_path) # DEBUG

        conflicting_sex_are_identical = compare_bgzed_txts(path_local(self.ref_conflicting_sex_path),
                                                           path_local(conflicting_sex_path))

        f_stat_outliers_are_identical = compare_bgzed_txts(path_local(self.ref_sex_annotation_f_stat_outliers_path),
                                                           path_local(sex_annotation_f_stat_outliers_path))


        self.assertTrue(conflicting_sex_are_identical and f_stat_outliers_are_identical)

    def test_2_2_1_prune_mt(self):
        # read reference output of the Step 2.1 impute_sex()
        ref_output_sex_mt = hl.read_matrix_table(self.ref_output_sex_mt_path)

        # run function to test
        pruned_mt = qc_step_2_2.prune_mt(ref_output_sex_mt, self.config)

        # compare output to reference
        output_mts_are_identical = compare_matrixtables(self.ref_mt_ldpruned_path, self.mt_ldpruned_path)

        self.assertTrue(output_mts_are_identical)

    def test_2_2_2_run_pc_relate(self):
        # read reference outputs of the Step 2.2 prune_mt()
        ref_mt_ldpruned = hl.read_matrix_table(self.ref_mt_ldpruned_path)

        # run function to test
        qc_step_2_2.run_pc_relate(ref_mt_ldpruned, self.config)

        # compare outputs to reference
        relatedness_ht_path = self.config['step2']['prune_pc_relate']['relatedness_ht_file']
        samples_to_remove_path = self.config['step2']['prune_pc_relate']['samples_to_remove_file']
        scores_path = self.config['step2']['prune_pc_relate']['scores_file']
        
        output_relatedness_ht_identical = compare_tables(relatedness_ht_path, self.ref_relatedness_ht_path)
        output_samples_to_remove_identical = compare_tables(samples_to_remove_path, self.ref_samples_to_remove_path)
        output_scores_identical = compare_tables(scores_path, self.ref_scores_path)

        self.assertTrue(output_relatedness_ht_identical and output_samples_to_remove_identical and output_scores_identical)

    def test_2_2_3_run_population_pca(self):
        # read reference outputs of Step 2.3
        # run function to test
        ref_mt_ldpruned = hl.read_matrix_table(self.ref_mt_ldpruned_path)
        ref_samples_to_remove = hl.read_table(self.ref_samples_to_remove_path)
        qc_step_2_2.run_population_pca(ref_mt_ldpruned, ref_samples_to_remove, self.config)
        
        # compare outputs to reference
        # pca_plots_identical = compare_txts(self.pca_output_path, self.ref_pca_output_path) # TODO: make this robust
        pca_scores_identical = compare_tables(self.mt_pca_scores_path, self.ref_mt_pca_scores_path)
        pca_loadings_identical = compare_tables(self.mt_pca_loadings_path, self.ref_mt_pca_loadings_path)
        pca_mt_identical = compare_matrixtables(self.pca_mt_path, self.ref_pca_mt_path)
        plink_identical = compare_plinks(self.plink_path + '.bed', self.plink_path + '.bim', self.plink_path + '.fam', 
                self.ref_plink_path + '.bed', self.ref_plink_path + '.bim', self.ref_plink_path + '.fam')

        self.assertTrue(pca_scores_identical and pca_loadings_identical and 
                        pca_mt_identical and plink_identical)

    def test_2_3_1_create_1kg_mt(self):
        # run function to test
        qc_step_2_3.create_1kg_mt(self.config)
        kg_wes_regions_identical = compare_matrixtables(self.kg_wes_regions, self.ref_kg_wes_regions)

        self.assertTrue(kg_wes_regions_identical)

    def test_2_3_2_merge_with_1kg(self):
        # use reference outputs of Step 2.2 prune_mt()
        ref_mt_ldpruned = hl.read_matrix_table(self.ref_mt_ldpruned_path)
        # run function to test
        qc_step_2_3.merge_with_1kg(ref_mt_ldpruned, self.config)
        merged_mt_identical = compare_matrixtables(self.merged_mt_path, self.ref_merged_mt_path)

        self.assertTrue(merged_mt_identical)

    def test_2_3_3_annotate_and_filter(self):
        # use reference outputs of Step 2.3 merge_with_1kg()
        ref_merged_mt = hl.read_matrix_table(self.ref_merged_mt_path)
        # run function to test
        qc_step_2_3.annotate_and_filter(ref_merged_mt, self.config)
        filtered_mt_identical = compare_matrixtables(self.filtered_mt_path, self.ref_filtered_mt_path)

        self.assertTrue(filtered_mt_identical)

    def test_2_3_4_run_pca(self):
        # use reference outputs of Step 2.3 annotate_and_filter()
        ref_filtered_mt = hl.read_matrix_table(self.ref_filtered_mt_path)
        # run function to test
        qc_step_2_3.run_pca(ref_filtered_mt, self.config)

        pca_scores_identical = compare_tables(self.pca_scores_path, self.ref_pca_scores_path)
        pca_loadings_identical = compare_tables(self.pca_loadings_path, self.ref_pca_loadings_path)
        pca_evals_identical = compare_txts(self.pca_evals_path, self.ref_pca_evals_path)

        self.assertTrue(pca_scores_identical and pca_loadings_identical and pca_evals_identical)

    def test_2_3_5_predict_pops(self):
        # use reference outputs of Step 2.3 run_pca()
        # run function to test
        qc_step_2_3.predict_pops(self.config)

        pop_ht_identical = compare_tables(self.pop_ht_file, self.ref_pop_ht_file)
        pop_ht_tsv_identical = compare_txts(path_local(self.pop_ht_tsv), path_local(self.ref_pop_ht_tsv))

        self.assertTrue(pop_ht_identical and pop_ht_tsv_identical)

    # vk11 : unprocessed after this line

    def test_2_4_1_annotate_mt(self):
        qc_step_2_4.annotate_mt(self.ref_mt_file, self.ref_pop_ht_file, self.annotated_mt_file)

        annotated_mts_identical = compare_matrixtables(self.annotated_mt_file, self.ref_annotated_mt_file)

        self.assertTrue(annotated_mts_identical)

    def test_2_4_2_stratified_sample_qc(self):
        qc_step_2_4.stratified_sample_qc(self.ref_annotated_mt_file, self.mt_qc_outfile,
                self.ht_qc_cols_outfile, self.qc_filter_file, self.config)

        mt_qc_identical = compare_matrixtables(self.mt_qc_outfile, self.ref_mt_qc_outfile)
        ht_qc_cols_identical = compare_tables(self.ht_qc_cols_outfile, self.ref_ht_qc_cols_outfile)
        qc_filter_identical = compare_tables(self.qc_filter_file, self.ref_qc_filter_file)
        out_text_identical = compare_bgzed_txts(path_local(self.output_text_file), 
                path_local(self.ref_output_text_file))

        unique_hail_id_replace = [[r'__uid_\d+\n', '']]
        out_globals_json_identical = compare_txts(path_local(self.output_globals_json), 
                path_local(self.ref_output_globals_json), replace_strings=unique_hail_id_replace)

        self.assertTrue(mt_qc_identical and ht_qc_cols_identical and 
                qc_filter_identical and out_text_identical and out_globals_json_identical)

    def test_2_5_1_remove_sample_qc_fails(self):
        qc_step_2_5.remove_sample_qc_fails(self.ref_annotated_mt_file, self.ref_qc_filter_file,
                self.samples_failing_qc_file, self.sample_qc_filtered_mt_file)

        sample_qc_filtered_identical = compare_matrixtables(self.sample_qc_filtered_mt_file, 
                self.ref_sample_qc_filtered_mt_file)
        samples_failing_qc_identical = compare_bgzed_txts(path_local(self.samples_failing_qc_file),
                path_local(self.ref_samples_failing_qc_file))

        self.assertTrue(sample_qc_filtered_identical and samples_failing_qc_identical)


    @classmethod
    def tearDownClass(cls):
        # clean up the directories created by the QC steps
        if hfs.is_dir(cls.mtdir):
            hfs.rmtree(cls.mtdir) 
        if hfs.is_dir(cls.annotdir):
            hfs.rmtree(cls.annotdir)
        if os.path.exists(cls.plotdir):
            sh.rmtree(cls.plotdir)
        super(TestQCSteps, cls).tearDownClass()


if __name__ == '__main__':
    unittest.main()

