import os
import re
import gzip
import unittest
import importlib

import hail as hl
import shutil as sh
import hailtop.fs as hfs

from typing import Union
from pyspark import SparkContext
from utils.utils import download_test_data_using_files_list
from utils.config import path_local, path_spark, _process_cvars_in_config


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
        # both fields not defined
        field_not_defined = (~hl.is_defined(struct1[field]) & ~hl.is_defined(struct2[field]))

        # both fields defined
        field_defined = hl.is_defined(struct1[field]) & hl.is_defined(struct2[field])

        # different checks for differenet field dtypes
        if struct1[field].dtype == hl.dtype('float'):
            # both fields not a number
            field_nan = hl.is_nan(struct1[field]) & hl.is_nan(struct2[field])
            expr = (field_defined & hl.approx_equal(struct1[field], struct2[field])) | field_not_defined | field_nan
        elif struct1[field].dtype == hl.dtype('array<float>'):
            # element_wise_expr = hl.map(lambda a, b: (hl.is_defined(a) & hl.is_defined(b) & hl.approx_equal(a, b)) | (~hl.is_defined(a) & ~hl.is_defined(b)), struct1[field], struct2[field])

            element_wise_expr = hl.map(lambda a, b: (hl.is_defined(a) &
                                                     hl.is_defined(b) &
                                                     hl.approx_equal(a, b)) | 
                                                     (~hl.is_defined(a) & 
                                                      ~hl.is_defined(b)) |
                                                     (hl.is_nan(a) &
                                                      hl.is_nan(b)),
                                       struct1[field], struct2[field])
            expr = hl.all(element_wise_expr)
        else:
            # TODO: implement comparison for other dtypes
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
    mt1 = hl.import_plink(bed=path_spark(bed1_path), bim=path_spark(bim1_path), fam=path_spark(fam1_path), reference_genome='GRCh38')
    mt2 = hl.import_plink(bed=path_spark(bed2_path), bim=path_spark(bim2_path), fam=path_spark(fam2_path), reference_genome='GRCh38')

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
        print('MatrixTables are equal')
        return True
    else:
        print(f'MatrixTables are not equal: {num_differences} differing entries found')
        return False

def compare_matrixtables(mt1_path: Union[str, hl.MatrixTable], mt2_path: Union[str, hl.MatrixTable]) -> bool:
    """
    Compare the contents of two Hail MatrixTables.

    Parameters:
    mt1_path (str): Path to the first MatrixTable.
    mt2_path (str): Path to the second MatrixTable.

    Returns:
    bool: True if the MatrixTables are equal, False otherwise.
    """
    # Load the MatrixTables
    if isinstance(mt1_path, str):
        mt1 = hl.read_matrix_table(path_spark(mt1_path))
    elif isinstance(mt1_path, hl.MatrixTable):
        mt1 = mt1_path
    else:
        raise TypeError(f"Expected 'str' or 'hl.Table', but got {type(mt1_path)}")
    
    if isinstance(mt2_path, str):
        mt2 = hl.read_matrix_table(path_spark(mt2_path))
    elif isinstance(mt2_path, hl.MatrixTable):
        mt2 = mt2_path
    else:
        raise TypeError(f"Expected 'str' or 'hl.Table', but got {type(mt2_path)}")

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
        print('MatrixTables are equal')
        return True
    else:
        print(f'MatrixTables {mt1_path}\n{mt2_path}\nare not equal: {num_differences} differing rows found')
        return False

def compare_entries_to_other_ht(this_row_value, row_key, other_ht):
    """
    Compare rows of two Tables.
    
    Returns:
    bool: True if entries are equal, False otherwise.
    """
    other_row_value = other_ht[row_key]
    return ~compare_structs(this_row_value, other_row_value)

def compare_tables(ht1_path: Union[str, hl.Table], ht2_path: Union[str, hl.Table]) -> bool:
    """
    Compare the contents of two Hail Tables.

    Parameters:
    ht1_path (str): Path to the first Table.
    ht2_path (str): Path to the second Table.

    Returns:
    bool: True if the Tables are equal, False otherwise.
    """
    # Load the Tables
    if isinstance(ht1_path, str):
        ht1 = hl.read_table(path_spark(ht1_path))
    elif isinstance(ht1_path, hl.Table):
        ht1 = ht1_path
    else:
        raise TypeError(f"Expected 'str' or 'hl.Table', but got {type(ht1_path)}")
    
    if isinstance(ht2_path, str):
        ht2 = hl.read_table(path_spark(ht2_path))
    elif isinstance(ht2_path, hl.Table):
        ht2 = ht2_path
    else:
        raise TypeError(f"Expected 'str' or 'hl.Table', but got {type(ht2_path)}")

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
        print(f'Tables {ht1_path}\n{ht2_path}\nare equal')
        return True
    else:
        print(f'Tables {ht1_path}\n{ht2_path}\nare not equal: {num_differences} differing rows found')
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
        print(f'Files {path_1}\n{path_2}\n are equal')
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
        print(f'Files {path_1}\n{path_2}\n are equal')
    else:
        print(f'Files {path_1}\n{path_2}\n are not equal')

    return contents_1 == contents_2

# DEBUG: might be useful https://spark.apache.org/docs/latest/api/python/getting_started/testing_pyspark.html
# TODO: clean up hail logs

# /path/to/wes-qc must be in PYTHONPATH
qc_step_1_1 = importlib.import_module("1-import_data.1-import_gatk_vcfs_to_hail")

qc_step_2_1 = importlib.import_module("2-sample_qc.1-hard_filters_sex_annotation")
qc_step_2_2 = importlib.import_module("2-sample_qc.2-prune_related_samples")
qc_step_2_3 = importlib.import_module("2-sample_qc.3-population_pca_prediction")
qc_step_2_4 = importlib.import_module("2-sample_qc.4-find_population_outliers")
qc_step_2_5 = importlib.import_module("2-sample_qc.5-filter_fail_sample_qc")

qc_step_3_1 = importlib.import_module("3-variant_qc.variant_qc_non_trios.1-generate_truth_sets_non_trios")
qc_step_3_2 = importlib.import_module("3-variant_qc.variant_qc_non_trios.2-create_rf_ht_non_trios")
qc_step_3_3 = importlib.import_module("3-variant_qc.3-train_rf")
qc_step_3_4 = importlib.import_module("3-variant_qc.4-apply_rf")
qc_step_3_5 = importlib.import_module("3-variant_qc.variant_qc_non_trios.5-annotate_ht_after_rf_no_trios")
qc_step_3_6 = importlib.import_module("3-variant_qc.variant_qc_non_trios.6-rank_and_bin_no_trios")
qc_step_3_7 = importlib.import_module("3-variant_qc.variant_qc_non_trios.7-plot_rf_output_no_trios")
qc_step_3_8 = importlib.import_module("3-variant_qc.8-select_thresholds")
qc_step_3_9 = importlib.import_module("3-variant_qc.9-filter_mt_after_variant_qc")

qc_step_4_1 = importlib.import_module("4-genotype_qc.1-apply_hard_filters")
qc_step_4_1a = importlib.import_module("4-genotype_qc.1a-apply_range_of_hard_filters")
qc_step_4_2 = importlib.import_module("4-genotype_qc.2-counts_per_sample")
qc_step_4_3 = importlib.import_module("4-genotype_qc.3-export_vcfs")
qc_step_4_3a = importlib.import_module("4-genotype_qc.3a-export_vcfs_range_of_hard_filters")
qc_step_4_3b = importlib.import_module("4-genotype_qc.3b-export_vcfs_stingent_filters")


class HailTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # outdir for the outputs of the regression tests
        cls.test_outdir_path = os.path.dirname(os.path.realpath(__file__))

        # tmp dir for Hail
        cls.tmp_dir = path_spark(os.path.join(cls.test_outdir_path, 'tmp_test'))

        cls.sc = SparkContext.getOrCreate()

        # initialise Hail in local mode to run the test even without cluster
        hl.init(sc=cls.sc, app_name="HailTest", master="local[*]",
                default_reference="GRCh38", tmp_dir=cls.tmp_dir)

    @classmethod
    def tearDownClass(cls):
        # clean up tmp_dir
        hfs.rmtree(cls.tmp_dir)
        hl.stop()
        cls.sc.stop()


# list with the test files must be located in the test directory in the repo
TEST_FILES_LIST = '../test_files_list_in_bucket.txt'

class TestQCSteps(HailTestCase):
    @classmethod
    def setUpClass(cls):
        # set up Hail from the base class
        super(TestQCSteps, cls).setUpClass()

        # `tests/` folder in the repo
        cls.test_suite_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

        cls.test_data_download_path = os.path.join(cls.test_suite_path, 'test_data_individual_files')

        # static test data used in multiple functions 
        cls.control_data_dir = os.path.join(cls.test_data_download_path, 'control_set_small')
        cls.resourcedir = os.path.join(cls.test_data_download_path, 'resources')
        cls.onekg_resourcedir = os.path.join(cls.resourcedir, 'mini_1000G')
        cls.training_sets = os.path.join(cls.test_data_download_path, 'training_sets')
        variant_qc_random_forest_path = os.path.join(cls.test_data_download_path, 'variant_qc_random_forest')

        # download test data from the bucket
        print(f'Downloading control set from the s3 bucket')
        download_test_data_using_files_list(TEST_FILES_LIST, cls.test_data_download_path)


        # define parameters needed for functions to be tested

        # TODO: separate input resources and annotations from output
        cls.ref_dataset_path = os.path.join(cls.test_data_download_path, 'unit_tests', 'reference_output_data')
        cls.ref_mtdir = os.path.join(cls.ref_dataset_path, 'matrixtables')
        cls.ref_resourcedir = os.path.join(cls.ref_dataset_path, 'resources')
        cls.ref_annotdir = os.path.join(cls.ref_dataset_path, 'annotations')
        cls.ref_plots_dir = os.path.join(cls.ref_dataset_path, 'plots_dir') # TODO: add to the bucket  


        # ===== QC General Params ===== #
        # output dirs for the regression tests outputs
        cls.out_mtdir = os.path.join(cls.test_outdir_path, 'matrixtables_test')
        cls.out_annotdir = os.path.join(cls.test_outdir_path, 'annotations_test')
        cls.out_resourcedir = os.path.join(cls.test_outdir_path, 'resources_test')
        cls.out_plots_dir = os.path.join(cls.test_outdir_path, 'plots_test')
        cls.out_var_qc_rf_dir = os.path.join(cls.test_suite_path, 'var_qc_rf_dir_test')

        # add general params into the config object
        config = dict()
        config['cvars'] = {
            'tmpdir': 'general.tmp_dir',
            'anndir': 'general.annotation_dir',
            'mtdir': 'general.matrixtables_dir',
            'resdir': 'general.resource_dir',
            '1kg_resdir': 'general.onekg_resource_dir',
            'pltdir': 'general.plots_dir',
            'traindir': 'general.training_sets_dir',
            'rfdir': 'general.var_qc_rf_dir'
        }

        config['general'] = dict(
            tmp_dir = cls.tmp_dir,
            annotation_dir = cls.out_annotdir,
            matrixtables_dir = cls.out_mtdir,
            resource_dir = cls.resourcedir,
            plots_dir = cls.out_plots_dir,
            onekg_resource_dir = cls.onekg_resourcedir, 
            training_sets_dir = cls.training_sets, 
            var_qc_rf_dir = cls.out_var_qc_rf_dir
        )

        # ===== QC Step 1.1 ===== #
        # inputs: `cls.control_data_dir`
        cls.vcf_header = '' # not available for the test dataset
        # outputs
        cls.output_mt_path = os.path.join(cls.out_mtdir, 'gatk_unprocessed.mt')
        # reference outputs
        cls.ref_mt_file = os.path.join(cls.ref_mtdir, 'gatk_unprocessed.mt')

        # add step 1 parameters into the config obj
        config['step1'] = dict(
            gatk_vcf_header_infile = cls.vcf_header,
            gatk_vcf_indir = cls.control_data_dir,
            gatk_mt_outfile = '{mtdir}/gatk_unprocessed.mt',
        )

        # ===== QC Step 2.1 ===== #
        # apply_hard_filters_outputs()
        # reference inputs
        cls.ref_unfiltered_mt_path = os.path.join(cls.ref_mtdir, 'gatk_unprocessed.mt')
        # outputs
        cls.output_filtered_mt_path = os.path.join(cls.out_mtdir, 'mt_hard_filters_annotated.mt')
        # reference outputs
        cls.ref_output_filtered_mt_path = os.path.join(cls.ref_mtdir, 'mt_hard_filters_annotated.mt')

        # impute_sex()
        # reference inputs: `cls.ref_output_filtered_mt_path`
        # outputs
        cls.output_sex_annotated_path = os.path.join(cls.out_annotdir, 'sex_annotated.sex_check.txt.bgz')
        cls.output_sex_mt_path = os.path.join(cls.out_mtdir, 'mt_sex_annotated.mt')
        # reference outputs
        cls.ref_output_sex_annotated_path = os.path.join(cls.ref_annotdir, 'sex_annotated.sex_check.txt.bgz')
        cls.ref_output_sex_mt_path = os.path.join(cls.ref_mtdir, 'mt_sex_annotated.mt')

        # identify_inconsistencies()
        # reference inputs: `cls.ref_output_sex_mt_path`
        # outputs
        cls.conflicting_sex_path = os.path.join(cls.out_annotdir, 'conflicting_sex.txt.bgz')
        cls.sex_annotation_f_stat_outliers_path = os.path.join(cls.out_annotdir, 'sex_annotation_f_stat_outliers.txt.bgz')
        # reference outputs
        cls.ref_conflicting_sex_path = os.path.join(cls.ref_annotdir, 'conflicting_sex.txt.bgz')
        cls.ref_sex_annotation_f_stat_outliers_path = os.path.join(cls.ref_annotdir, 'sex_annotation_f_stat_outliers.txt.bgz')

        # add step 2.1 parameters into the config obj
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
        cls.mt_ldpruned_path = os.path.join(cls.out_mtdir, 'mt_ldpruned.mt')
        # reference outputs
        cls.ref_mt_ldpruned_path = os.path.join(cls.ref_mtdir, 'mt_ldpruned.mt')

        # run_pc_relate()
        # reference inputs: `cls.ref_mt_ldpruned_path`
        # outputs
        cls.relatedness_ht_path = os.path.join(cls.out_mtdir, 'mt_relatedness.ht')
        cls._path_path = os.path.join(cls.out_mtdir, 'mt_related_samples_to_remove.ht')
        cls.scores_path = os.path.join(cls.out_mtdir, 'mt_pruned.pca_scores.ht')
        # reference outputs
        cls.ref_relatedness_ht_path = os.path.join(cls.ref_mtdir, 'mt_relatedness.ht')
        cls.ref_samples_to_remove_path = os.path.join(cls.ref_mtdir, 'mt_related_samples_to_remove.ht')
        cls.ref_scores_path = os.path.join(cls.ref_mtdir, 'mt_pruned.pca_scores.ht')

        # run_population_pca()
        # reference inputs: `cls.ref_mt_ldpruned_path`, `cls.ref_samples_to_remove_path`
        # outputs
        cls.plotdir = os.path.join(cls.test_outdir_path, 'plots_test')
        cls.plink_path = os.path.join(cls.out_mtdir, 'mt_unrelated.plink')
        cls.mt_pca_scores_path = os.path.join(cls.out_mtdir, 'mt_pca_scores.ht')
        cls.mt_pca_loadings_path = os.path.join(cls.out_mtdir, 'mt_pca_loadings.ht')
        cls.pca_output_path = os.path.join(cls.plotdir, 'pca.html')
        cls.pca_mt_path = os.path.join(cls.out_mtdir, 'mt_pca.mt')
        # reference outputs
        cls.ref_plink_path = os.path.join(cls.ref_mtdir, 'mt_unrelated.plink')
        cls.ref_mt_pca_scores_path = os.path.join(cls.ref_mtdir, 'mt_pca_scores.ht')
        cls.ref_mt_pca_loadings_path = os.path.join(cls.ref_mtdir, 'mt_pca_loadings.ht')
        cls.ref_pca_output_path = path_local(os.path.join(cls.ref_dataset_path, 'plots', 'pca.html')) # DEBUG: not used
        cls.ref_pca_mt_path = os.path.join(cls.ref_mtdir, 'mt_pca.mt')

        config['step2']['prune'] = {
            'pruned_mt_file': '{mtdir}/mt_ldpruned.mt',
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
        cls.kg_wes_regions = os.path.join(cls.out_mtdir, 'kg_wes_regions.mt')
        # reference outputs
        cls.ref_kg_wes_regions = os.path.join(cls.ref_mtdir, 'kg_wes_regions.mt')

        # merge_with_1kg()
        # reference inputs: `cls.ref_mt_ldpruned_path`
        # outputs
        cls.merged_mt_path = os.path.join(cls.out_mtdir, 'merged_with_1kg.mt')
        # reference outputs
        cls.ref_merged_mt_path = os.path.join(cls.ref_mtdir, "merged_with_1kg.mt")

        # annotate_and_filter()
        # reference inputs: `cls.ref_merged_mt_path`, `cls.resourcedir`
        # outputs
        cls.filtered_mt_path = os.path.join(cls.out_mtdir, 'merged_with_1kg_filtered.mt')
        # reference outputs
        cls.ref_filtered_mt_path = os.path.join(cls.ref_mtdir, "merged_with_1kg_filtered.mt")

        # run_pca()
        # reference inputs: `cls.ref_filtered_mt_path`
        # outputs
        cls.pca_scores_path = os.path.join(cls.out_mtdir, 'pca_scores_after_pruning.ht')
        cls.pca_loadings_path = os.path.join(cls.out_mtdir, 'pca_loadings_after_pruning.ht')
        cls.pca_evals_path = path_local(os.path.join(cls.out_mtdir, 'pca_evals_after_pruning.txt'))
        # reference outputs
        cls.ref_pca_scores_path = os.path.join(cls.ref_mtdir, 'pca_scores_after_pruning.ht')
        cls.ref_pca_loadings_path = os.path.join(cls.ref_mtdir, 'pca_loadings_after_pruning.ht')
        cls.ref_pca_evals_path = path_local(os.path.join(cls.ref_mtdir, 'pca_evals_after_pruning.txt'))

        # predict_pops()
        # reference inputs: `cls.ref_pca_scores_path`
        # outputs
        cls.pop_ht_file = os.path.join(cls.out_mtdir, 'pop_assignments.ht')
        cls.pop_ht_tsv = os.path.join(cls.out_mtdir, 'pop_assignments.tsv')
        # reference outputs
        cls.ref_pop_ht_file = os.path.join(cls.ref_mtdir, 'pop_assignments.ht')
        cls.ref_pop_ht_tsv = os.path.join(cls.ref_mtdir, 'pop_assignments.tsv')

        config['step2']['create_1kg_mt'] = {
            'indir': '{1kg_resdir}',
            'vcfheader': '{1kg_resdir}/header_20201028.txt',
            'mt_out_file' : '{mtdir}/kg_wes_regions.mt'
        }

        config['step2']['merge_with_1_kg'] = {
            'kg_mt_file' : cls.ref_kg_wes_regions,
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
            'pca_scores_file' : cls.ref_pca_scores_path,
            'gnomad_pc_n_estimators' : 100,
            'gnomad_prop_train' : 0.8,
            'gnomad_min_prob' : 0.5,
            'pop_ht_outfile' : '{mtdir}/pop_assignments.ht',
            'pop_ht_outtsv' : '{mtdir}/pop_assignments.tsv'}

        # ===== QC Step 2.4 ===== #
        # annotate_mt()
        # reference inputs: `cls.ref_pop_ht_file`, `cls.ref_pop_ht_tsv`, `cls.ref_mt_file`
        # outputs
        cls.annotated_mt_file = os.path.join(cls.out_mtdir, 'gatk_unprocessed_with_pop.mt')
        # reference outputs
        cls.ref_annotated_mt_file = os.path.join(cls.ref_mtdir, 'gatk_unprocessed_with_pop.mt')

        # stratified_sample_qc()
        # reference inputs: `cls.ref_annotated_mt_file`, `cls.ref_annotdir`
        # outputs
        cls.mt_qc_outfile = os.path.join(cls.out_mtdir, 'mt_pops_sampleqc.mt')
        cls.ht_qc_cols_outfile = os.path.join(cls.out_mtdir, 'mt_pops_sampleqc.ht')
        cls.qc_filter_file = os.path.join(cls.out_mtdir, 'mt_pops_QC_filters.ht')
        cls.output_text_file = os.path.join(cls.out_annotdir, 'sample_qc_by_pop.tsv.bgz')
        cls.output_globals_json = os.path.join(cls.out_annotdir, 'sample_qc_by_pop.globals.json')

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
        cls.sample_qc_filtered_mt_file = os.path.join(cls.out_mtdir, 'mt_pops_QC_filters_after_sample_qc.mt')
        cls.samples_failing_qc_file = os.path.join(cls.out_annotdir, 'samples_failing_qc.tsv.bgz')
        # reference outputs
        cls.ref_sample_qc_filtered_mt_file = os.path.join(cls.ref_mtdir, 'mt_pops_QC_filters_after_sample_qc.mt')
        cls.ref_samples_failing_qc_file = os.path.join(cls.ref_annotdir, 'samples_failing_qc.tsv.bgz')

        config['step2']['remove_sample_qc_fails'] = {
            'samples_failing_qc_file' : '{anndir}/samples_failing_qc.tsv.bgz',
            'sample_qc_filtered_mt_file' : '{mtdir}/mt_pops_QC_filters_after_sample_qc.mt'
        }

        # ===== QC Step 3.1 non-trios ===== #
        # get_truth_ht()
        # reference inputs:
        # resource inputs:
        cls.omni = os.path.join(cls.training_sets, '1000G_omni2.5.hg38.ht')
        cls.mills = os.path.join(cls.training_sets, 'Mills_and_1000G_gold_standard.indels.hg38.ht')
        cls.thousand_genomes = os.path.join(cls.training_sets, '1000G_phase1.snps.high_confidence.hg38.ht')
        cls.hapmap = os.path.join(cls.training_sets, 'hapmap_3.3.hg38.ht')
        # outputs
        cls.truth_ht = os.path.join(cls.out_resourcedir, 'truthset_table.ht')
        # reference outputs
        cls.ref_truth_ht = os.path.join(cls.ref_resourcedir, 'truthset_table.ht')

        # split_multi_and_var_qc()
        # reference inputs: `cls.ref_sample_qc_filtered_mt_file`
        # outputs
        cls.mt_varqc = os.path.join(cls.out_mtdir, 'mt_varqc.mt')
        cls.mt_varqc_splitmulti = os.path.join(cls.out_mtdir, 'mt_varqc_splitmulti.mt') # TODO: implement comparison of any structures in Matrix Tables
        # reference outputs
        cls.ref_mt_varqc = os.path.join(cls.ref_mtdir, 'mt_varqc.mt')
        cls.ref_mt_varqc_splitmulti = os.path.join(cls.ref_mtdir, 'mt_varqc_splitmulti.mt') # TODO: implement comparison of any structures in Matrix Tables

        # create_inbreeding_ht_with_ac_and_allele_data()
        # reference inputs: `cls.ref_mt_varqc`
        # outputs
        cls.inbreeding_htoutfile = os.path.join(cls.out_mtdir, 'inbreeding.ht') # TODO: implement comparison of any structures in Hail Tables
        cls.qc_ac_htoutfile = os.path.join(cls.out_mtdir, 'qc_ac.ht')
        cls.allele_data_htoutfile = os.path.join(cls.out_mtdir, 'allele_data.ht')
        # reference outputs
        cls.ref_inbreeding_htoutfile = os.path.join(cls.ref_mtdir, 'inbreeding.ht') # TODO: implement comparison of any structures in Hail Tables
        cls.ref_qc_ac_htoutfile = os.path.join(cls.ref_mtdir, 'qc_ac.ht')
        cls.ref_allele_data_htoutfile = os.path.join(cls.ref_mtdir, 'allele_data.ht')

        # ===== QC Step 3.2 non-trios ===== #
        # create_rf_ht()
        # reference inputs: `cls.ref_mt_varqc_splitmulti`, ``
        # outputs
        cls.htoutfile_rf_all_cols = os.path.join(cls.out_mtdir, 'ht_for_RF_all_cols.ht')
        cls.htoutfile_rf_var_type_all_cols = os.path.join(cls.out_mtdir, 'ht_for_RF_by_variant_type_all_cols.ht')
        # reference outputs
        cls.ref_htoutfile_rf_all_cols = os.path.join(cls.ref_mtdir, 'ht_for_RF_all_cols.ht')
        cls.ref_htoutfile_rf_var_type_all_cols = os.path.join(cls.ref_mtdir, 'ht_for_RF_by_variant_type_all_cols.ht')

        config['step3'] = dict()
        config['step3']['create_rf_ht'] = {
            'fail_hard_filters_QD_less_than': 2,
            'fail_hard_filters_FS_greater_than': 60,
            'fail_hard_filters_MQ_less_than': 30,
        }


       
        # ===== QC Step 4.1 ===== #
        # filter_mt()
        # reference inputs: `cls.ref_mt_after_var_qc`
        cls.ref_mt_after_var_qc = os.path.join(cls.ref_mtdir, 'mt_after_var_qc.mt')
        # outputs
        cls.mtfile_filtered = os.path.join(cls.out_mtdir, 'mt_after_var_qc_hard_filter_gt.mt')
        # reference outputs
        cls.ref_mtfile_filtered = os.path.join(cls.ref_mtdir, 'mt_after_var_qc_hard_filter_gt.mt')


        # ===== QC Step 4.2 ===== #
        # annotate_gnomad() # 
        # # reference inputs:
        # cls.ref_mt_hard_filter_combinations = os.path.join(cls.ref_mtdir, 'mt_hard_filter_combinations.mt')
        # # resource inputs:
        # cls.gnomad_htfile = os.path.join(cls.resourcedir, 'gnomad.exomes.r2.1.1.sites.liftover_grch38.ht')
        # # outputs
        # cls.annotated_with_gnomad = os.path.join(cls.out_mtdir, 'annotated_with_gnomad.mt') # not used as the is returned from the function, not written on disk
        # # reference outputs
        # cls.ref_annotated_with_gnomad = os.path.join(cls.ref_mtdir, 'annotated_with_gnomad.mt')

        # get_counts_per_cq()
        # reference inputs: `cls.ref_annotated_with_gnomad`
        # outputs
        cls.cqfile = os.path.join(cls.out_plots_dir, 'variant_counts_per_cq_post_qc.txt')
        # reference outputs
        cls.ref_cqfile = os.path.join(cls.ref_plots_dir, 'variant_counts_per_cq_post_qc.txt')

        # get_ca_fractions()
        # reference inputs: `cls.ref_annotated_with_gnomad`
        # outputs
        cls.cafile = os.path.join(cls.out_plots_dir, 'frac_ca_per_sample_post_qc_snv.txt')
        # reference outputs
        cls.ref_cafile = os.path.join(cls.ref_plots_dir, 'frac_ca_per_sample_post_qc_snv.txt')

        cls.config = _process_cvars_in_config(config)
        
    # QC Step 1.1
    def test_1_1_load_vcfs_to_mt(self):
        # run function to test
        qc_step_1_1.load_vcfs_to_mt(self.config)
        
        # compare the output to reference
        output_mt_path = self.config['step1']['gatk_mt_outfile']
        self.assertEqual(output_mt_path, self.output_mt_path)

        outputs_are_identical = compare_matrixtables(self.ref_mt_file, output_mt_path)
        self.assertTrue(outputs_are_identical)

    # QC Step 2.1 apply_hard_filters
    def test_2_1_1_apply_hard_filters(self):
        # read reference output of the step 1.1
        ref_mt_unfiltered = hl.read_matrix_table(path_spark(self.ref_unfiltered_mt_path))

        # run function to test
        mt_filtered = qc_step_2_1.apply_hard_filters(ref_mt_unfiltered, self.config)
        
        # compare the output to reference
        output_filtered_mt_path = self.config['step2']['sex_annotation_hard_filters']['filtered_mt_outfile']
        self.assertEqual(output_filtered_mt_path, self.output_filtered_mt_path)
        
        outputs_are_identical = compare_matrixtables(self.ref_output_filtered_mt_path, output_filtered_mt_path)
        self.assertTrue(outputs_are_identical)

    # QC Step 2.1 impute sex
    def test_2_1_2_impute_sex(self):
        # read reference output of the Step 2.1 apply_hard_filters()
        ref_mt_filtered = hl.read_matrix_table(path_spark(self.ref_output_filtered_mt_path))

        # run function to test
        sex_mt = qc_step_2_1.impute_sex(ref_mt_filtered, self.config)

        # compare the outputs to reference
        output_sex_annotated_path = self.config['step2']['impute_sex']['sex_ht_outfile']
        output_sex_mt_path = self.config['step2']['impute_sex']['sex_mt_outfile']
        self.assertEqual(output_sex_annotated_path, self.output_sex_annotated_path)
        self.assertEqual(output_sex_mt_path, self.output_sex_mt_path)

        self.assertEqual(path_local(output_sex_annotated_path), path_local(output_sex_annotated_path))

        output_mts_are_identical = compare_matrixtables(self.ref_output_sex_mt_path, output_sex_mt_path)
        output_txts_are_identical = compare_bgzed_txts(path_local(self.ref_output_sex_annotated_path),
                                                       path_local(output_sex_annotated_path))

        self.assertTrue(output_mts_are_identical and output_txts_are_identical)

    def test_2_1_3_identify_inconsistencies(self):
        # read reference output of the Step 2.1 impute_sex()
        ref_output_sex_mt = hl.read_matrix_table(path_spark(self.ref_output_sex_mt_path))

        # run function to test
        qc_step_2_1.identify_inconsistencies(ref_output_sex_mt, self.config)

        # compare outputs to reference
        conflicting_sex_path = self.config['step2']['sex_inconsistencies']['conflicting_sex_report_file']
        sex_annotation_f_stat_outliers_path = self.config['step2']['sex_inconsistencies']['fstat_outliers_report_file']

        self.assertEqual(conflicting_sex_path, self.conflicting_sex_path)
        self.assertEqual(sex_annotation_f_stat_outliers_path, self.sex_annotation_f_stat_outliers_path)

        conflicting_sex_are_identical = compare_bgzed_txts(path_local(self.ref_conflicting_sex_path),
                                                           path_local(conflicting_sex_path))

        f_stat_outliers_are_identical = compare_bgzed_txts(path_local(self.ref_sex_annotation_f_stat_outliers_path),
                                                           path_local(sex_annotation_f_stat_outliers_path))


        self.assertTrue(conflicting_sex_are_identical and f_stat_outliers_are_identical)

    def test_2_2_1_prune_mt(self):
        # read reference output of the Step 2.1 impute_sex()
        ref_output_sex_mt = hl.read_matrix_table(path_spark(self.ref_output_sex_mt_path))

        # run function to test
        pruned_mt = qc_step_2_2.prune_mt(ref_output_sex_mt, self.config)

        # compare output to reference
        output_mts_are_identical = compare_matrixtables(self.ref_mt_ldpruned_path, self.mt_ldpruned_path)

        self.assertTrue(output_mts_are_identical)

    def test_2_2_2_run_pc_relate(self):
        # read reference outputs of the Step 2.2 prune_mt()
        ref_mt_ldpruned = hl.read_matrix_table(path_spark(self.ref_mt_ldpruned_path))

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
        ref_mt_ldpruned = hl.read_matrix_table(path_spark(self.ref_mt_ldpruned_path))
        ref_samples_to_remove = hl.read_table(path_spark(self.ref_samples_to_remove_path))
        qc_step_2_2.run_population_pca(ref_mt_ldpruned, ref_samples_to_remove, self.config)
        
        # compare outputs to reference
        # pca_plots_identical = compare_txts(self.pca_output_path, self.ref_pca_output_path) # TODO: implement robust html comparison
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
        ref_mt_ldpruned = hl.read_matrix_table(path_spark(self.ref_mt_ldpruned_path))
        # run function to test
        qc_step_2_3.merge_with_1kg(ref_mt_ldpruned, self.config)
        merged_mt_identical = compare_matrixtables(self.merged_mt_path, self.ref_merged_mt_path)

        self.assertTrue(merged_mt_identical)

    def test_2_3_3_annotate_and_filter(self):
        # use reference outputs of Step 2.3 merge_with_1kg()
        ref_merged_mt = hl.read_matrix_table(path_spark(self.ref_merged_mt_path))
        # run function to test
        qc_step_2_3.annotate_and_filter(ref_merged_mt, self.config)
        filtered_mt_identical = compare_matrixtables(self.filtered_mt_path, self.ref_filtered_mt_path)

        self.assertTrue(filtered_mt_identical)

    def test_2_3_4_run_pca(self):
        # Use reference output of Step 2.3 annotate_and_filter() as input
        ref_filtered_mt = hl.read_matrix_table(path_spark(self.ref_filtered_mt_path))
        # Run function to test
        qc_step_2_3.run_pca(ref_filtered_mt, self.config)

        # Compare output files to the reference
        pca_scores_identical = compare_tables(self.pca_scores_path, self.ref_pca_scores_path)
        pca_loadings_identical = compare_tables(self.pca_loadings_path, self.ref_pca_loadings_path)
        pca_evals_identical = compare_txts(self.pca_evals_path, self.ref_pca_evals_path)

        # Test passes when all files match the references
        self.assertTrue(pca_scores_identical and pca_loadings_identical and pca_evals_identical)

    def test_2_3_5_predict_pops(self):
        # use reference outputs of Step 2.3 run_pca()
        # run function to test
        qc_step_2_3.predict_pops(self.config)

        pop_ht_identical = compare_tables(self.pop_ht_file, self.ref_pop_ht_file)
        pop_ht_tsv_identical = compare_txts(path_local(self.pop_ht_tsv), path_local(self.ref_pop_ht_tsv))

        self.assertTrue(pop_ht_identical and pop_ht_tsv_identical)

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

        unique_hail_id_replace = [[r'__uid_\d+\n', '']] # compare files up to a unique Hail file id
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


    # tests for 4-genotype_qc
    def test_4_1_1_filter_mt(self):
        # NOTE: when generating new reference make sure to use these arguments
        qc_step_4_1.filter_mt(self.ref_mt_after_var_qc, dp=5, gq=10, ab=0.2, 
                              mtfile_filtered=self.mtfile_filtered)

        mtfile_filetered_identical = compare_matrixtables(self.mtfile_filtered, self.ref_mtfile_filtered)

        self.assertTrue(mtfile_filetered_identical)

    # TODO: implement tests for the remaining functions listed below

    # test data has no trios, so only non-trios scripts are tested
    def test_3_non_trios_1_1_get_truth_ht(self):
        # TODO: generate reference
        qc_step_3_1.get_truth_ht(self.omni, self.mills, self.thousand_genomes,
                                 self.hapmap, self.truth_ht)
        

        truth_ht_identical = compare_tables(self.truth_ht, self.ref_truth_ht)

        self.assertTrue(truth_ht_identical)

    def test_3_non_trios_1_2_split_multi_and_var_qc(self):
        # TODO: generate reference
        qc_step_3_1.split_multi_and_var_qc(self.ref_sample_qc_filtered_mt_file, 
                                           self.mt_varqc, self.mt_varqc_splitmulti)
        
        mt_varqc_identical = compare_matrixtables(self.mt_varqc, self.ref_mt_varqc)
        # mt_varqc_splitmulti_identical = compare_matrixtables(self.mt_varqc_splitmulti, self.ref_mt_varqc_splitmulti) # TODO: implement comparison of any structures in Matrix Tables

        self.assertTrue(mt_varqc_identical)

    def test_3_non_trios_1_3_create_inbreeding_ht_with_ac_and_allele_data(self):
        qc_step_3_1.create_inbreeding_ht_with_ac_and_allele_data(self.ref_mt_varqc, self.inbreeding_htoutfile, 
                                                                 self.qc_ac_htoutfile, self.allele_data_htoutfile)
        
        # inbreeding_htoutfile_identical = compare_tables(self.inbreeding_htoutfile, self.ref_inbreeding_htoutfile) # TODO: implement comparison of any structures in Hail Tables
        qc_ac_htoutfile_identical = compare_tables(self.qc_ac_htoutfile, self.ref_qc_ac_htoutfile)
        allele_data_htoutfile_identical = compare_tables(self.allele_data_htoutfile, self.ref_allele_data_htoutfile)

        self.assertTrue(qc_ac_htoutfile_identical and allele_data_htoutfile_identical)

    def test_3_non_trios_2_1_create_rf_ht(self):
        qc_step_3_2.create_rf_ht(self.ref_mt_varqc_splitmulti, self.ref_truth_ht, 
                 self.ref_allele_data_htoutfile, self.ref_qc_ac_htoutfile, 
                 self.ref_inbreeding_htoutfile, self.htoutfile_rf_all_cols, 
                 self.htoutfile_rf_var_type_all_cols, self.config)
        
        htoutfile_rf_all_cols_identical = compare_tables(self.htoutfile_rf_all_cols, self.ref_htoutfile_rf_all_cols)
        htoutfile_rf_all_cols_var_type_all_colls_identical = compare_tables(self.htoutfile_rf_var_type_all_cols, self.ref_htoutfile_rf_var_type_all_cols)

        self.assertTrue(htoutfile_rf_all_cols_identical, htoutfile_rf_all_cols_var_type_all_colls_identical)

    @classmethod
    def tearDownClass(cls):
        # clean up the directories created by the QC steps
#        if hfs.is_dir(cls.out_mtdir):
#            hfs.rmtree(cls.out_mtdir) 
#        if hfs.is_dir(cls.out_annotdir):
#            hfs.rmtree(cls.out_annotdir)
#        if os.path.exists(cls.plotdir):
#            sh.rmtree(cls.plotdir)
        super(TestQCSteps, cls).tearDownClass()


if __name__ == '__main__':
    unittest.main()

