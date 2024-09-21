import os
import argparse 
from unittest.mock import patch 
import re
import unittest
import importlib
from typing import Optional
import hail as hl
import hailtop.fs as hfs
from pyspark import SparkContext
import subprocess
from utils.utils import download_test_data_from_s3

# for test config rendering
INTEGRATION_TESTS_DIR = '{INTEGRATION_TESTS_DIR}'
TEST_DATA_DIR = '{TEST_DATA_DIR}'
RESOURCES_DIR = '{RESOURCES_DIR}'
TRAINING_SETS_DIR = '{TRAINING_SETS_DIR}'
VARIANT_QC_RANDOM_FOREST_DIR = '{VARIANT_QC_RANDOM_FOREST_DIR}'

def render_config(path_to_template: str, test_data_dir: Optional[str],
                  resources_dir: Optional[str], training_sets_dir: Optional[str], variant_qc_random_forest_dir: Optional[str],
                  savefile: str='inputs_test_rendered.yaml'):
    """
    Read the config template and fill in the paths.
    """
    integration_tests_dir = os.path.dirname(os.path.abspath(__file__))

    # TODO agree on test and resources folder naming convention
    test_data_dir = test_data_dir if test_data_dir else os.path.join(integration_tests_dir, 'control_set') 
    resources_dir = resources_dir if resources_dir else os.path.join(integration_tests_dir, 'resources')
    training_sets_dir = training_sets_dir if training_sets_dir else os.path.join(integration_tests_dir, 'training_sets')
    variant_qc_random_forest_dir = variant_qc_random_forest_dir if variant_qc_random_forest_dir else os.path.join(integration_tests_dir, 'variant_qc_random_forest_dir')

    with open(path_to_template, 'r') as f:
        template = f.read()

    # TODO: make versatile
    template = template.replace(INTEGRATION_TESTS_DIR, integration_tests_dir)
    template = template.replace(TEST_DATA_DIR, test_data_dir)
    template = template.replace(RESOURCES_DIR, resources_dir)
    template = template.replace(TRAINING_SETS_DIR, training_sets_dir)
    template = template.replace(VARIANT_QC_RANDOM_FOREST_DIR, variant_qc_random_forest_dir)

    with open(savefile, 'w') as f:
        f.write(template)

# keep an eye on test execution order
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
qc_step_3_1 = importlib.import_module("3-variant_qc.1-generate_truth_sets")
qc_step_3_2 = importlib.import_module("3-variant_qc.2-create_rf_ht")
qc_step_3_3 = importlib.import_module("3-variant_qc.3-train_rf")
qc_step_3_4 = importlib.import_module("3-variant_qc.4-apply_rf")
qc_step_3_5 = importlib.import_module("3-variant_qc.5-annotate_ht_after_rf")
qc_step_3_6 = importlib.import_module("3-variant_qc.6-rank_and_bin")
qc_step_3_7 = importlib.import_module("3-variant_qc.7-plot_rf_output")
qc_step_3_8 = importlib.import_module("3-variant_qc.8-select_thresholds")
qc_step_3_9 = importlib.import_module("3-variant_qc.9-filter_mt_after_variant_qc")

# TEST_DATA_DOWNLOAD_URL = 'https://wes-qc-data.cog.sanger.ac.uk/all_test_data/test_data.zip' # moved to utils

class HailTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # get test suite path in the repo
        test_suite_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

        # specify paths to the test data and resources
        test_data_path = os.path.join(test_suite_path, 'control_set_small')
        resources_path = os.path.join(test_suite_path, 'resources')
        training_sets_path = os.path.join(test_suite_path, 'training_sets')
        variant_qc_random_forest_path = os.path.join(test_suite_path, 'variant_qc_random_forest')
        ref_data_path = os.path.join(test_suite_path, 'unit_tests')

        unzipped_path = os.path.join(test_suite_path, 'unzipped_data')
        unzipped_control_set_path = os.path.join(unzipped_path, 'control_set_small')
        unzipped_resources_path = os.path.join(unzipped_path, 'resources')
        unzipped_training_sets_path = os.path.join(unzipped_path, 'training_sets')

        # not used in integration tests, but prolly better download it as well
        unzipped_ref_data_path = os.path.join(unzipped_path, 'unit_tests', 'reference_output_data')


        test_data_dirs_to_move = {
            unzipped_control_set_path: test_suite_path,
            unzipped_resources_path: test_suite_path,
            unzipped_training_sets_path: test_suite_path,
            unzipped_ref_data_path: ref_data_path
        }

        download_test_data_from_s3(unzipped_path, test_data_dirs_to_move)

        smoke_test_dir_path = os.path.dirname(os.path.realpath(__file__))
        rendered_config_savefile = os.path.join(smoke_test_dir_path, 'integration_config_rendered.yaml')
        
        # # render test config from the template
        # render_config('inputs_test_template.yaml', test_data_path, resources_path) # TODO: make configurable
        render_config('new_config_test_template.yaml', test_data_path, resources_path, training_sets_path, variant_qc_random_forest_path,
                      savefile=rendered_config_savefile)

        # set up path to test config
        os.environ['WES_CONFIG'] = rendered_config_savefile

    @classmethod
    def tearDownClass(cls):
        # TODO: clean up logs
        pass

class IntegrationTests(HailTestCase):
    def test_1_1_import_data(self):
        try:
            qc_step_1_1.main()
        except Exception as e:
            self.fail(f'Step 1.1 failed with an exception: {e}')

    def test_2_1_sample_qc(self):
        try:
            qc_step_2_1.main()
        except Exception as e:
            self.fail(f'Step 2.1 failed with an exception: {e}')

    def test_2_2_sample_qc(self):
        try:
            qc_step_2_2.main()
        except Exception as e:
            self.fail(f'Step 2.2 failed with an exception: {e}')

    # mock cli arguments
    @patch('argparse.ArgumentParser.parse_args',
    return_value=argparse.Namespace(kg_to_mt=True, run=True, merge=True, filter=True, pca=True, assign_pops=True))
    def test_2_3_sample_qc(self, mock_args):
        try:
            qc_step_2_3.main()
        except Exception as e:
            self.fail(f'Step 2.3 failed with an exception: {e}')

    def test_2_4_sample_qc(self):
        try:
            qc_step_2_4.main()
        except Exception as e:
            self.fail(f'Step 2.4 failed with an exception: {e}')

    def test_2_5_sample_qc(self):
        try:
            qc_step_2_5.main()
        except Exception as e:
            self.fail(f'Step 2.5 failed with an exception: {e}')
    
    # mock cli arguments
    @patch('argparse.ArgumentParser.parse_args',
    return_value=argparse.Namespace(all=True))
    def test_3_1_variant_qc(self, mock_args):
        try:
            qc_step_3_1.main()
        except Exception as e:
            self.fail(f'Step 3.1 failed with an exception: {e}')
    
    def test_3_2_variant_qc(self):
        try:
            qc_step_3_2.main()
        except Exception as e:
            self.fail(f'Step 3.2 failed with an exception: {e}')

    def test_3_3_variant_qc(self):
        try:
            qc_step_3_3.main()
        except Exception as e:
            self.fail(f'Step 3.3 failed with an exception: {e}')

    # mock cli arguments
    @patch('argparse.ArgumentParser.parse_args',
    return_value=argparse.Namespace(runhash='beep')) # vk11: do not know what to set in here
    def test_3_4_variant_qc(self, mock_args):
        try:
            qc_step_3_4.main()
        except Exception as e:
            self.fail(f'Step 3.4 failed with an exception: {e}')

    # mock cli arguments
    @patch('argparse.ArgumentParser.parse_args',
    return_value=argparse.Namespace(runhash='beep')) 
    def test_3_5_variant_qc(self, mock_args):
        try:
            qc_step_3_5.main()
        except Exception as e:
            self.fail(f'Step 3.5 failed with an exception: {e}')

    # mock cli arguments
    @patch('argparse.ArgumentParser.parse_args',
    return_value=argparse.Namespace(runhash='beep')) 
    def test_3_6_variant_qc(self, mock_args):
        try:
            qc_step_3_6.main()
        except Exception as e:
            self.fail(f'Step 3.6 failed with an exception: {e}')

    # mock cli arguments
    @patch('argparse.ArgumentParser.parse_args',
    return_value=argparse.Namespace(runhash='beep')) 
    def test_3_7_variant_qc(self, mock_args):
        try:
            qc_step_3_7.main()
        except Exception as e:
            self.fail(f'Step 3.7 failed with an exception: {e}')

    # mock cli arguments
    @patch('argparse.ArgumentParser.parse_args',
    return_value=argparse.Namespace(runhash='beep', snv=0.1, indel=0.1)) # vk11: again, I do not know what values to set 
    def test_3_8_variant_qc(self, mock_args):
        try:
            qc_step_3_8.main()
        except Exception as e:
            self.fail(f'Step 3.8 failed with an exception: {e}')

    # mock cli arguments
    @patch('argparse.ArgumentParser.parse_args',
    return_value=argparse.Namespace(runhash='beep', snv=0.1, indel=0.1)) 
    def test_3_9_variant_qc(self):
        try:
            qc_step_3_9.main()
        except Exception as e:
            self.fail(f'Step 3.9 failed with an exception: {e}')

if __name__ == '__main__':
    unittest.main()

