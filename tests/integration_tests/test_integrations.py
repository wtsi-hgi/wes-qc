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

# for test config rendering
INTEGRATION_TESTS_DIR = '{INTEGRATION_TESTS_DIR}'
TEST_DATA_DIR = '{TEST_DATA_DIR}'
RESOURCES_DIR = '{RESOURCES_DIR}'

def render_config(path_to_template: str, test_data_dir: Optional[str],
                  resources_dir: Optional[str], savefile: str='inputs_test_rendered.yaml'):
    """
    Read the config template and fill in the paths.
    """
    integration_tests_dir = os.path.dirname(os.path.abspath(__file__))

    # TODO agree on test and resources folder naming convention
    test_data_dir = test_data_dir if test_data_dir else os.path.join(integration_tests_dir, 'control_set') 
    resources_dir = resources_dir if resources_dir else os.path.join(integration_tests_dir, 'resources')

    with open(path_to_template, 'r') as f:
        template = f.read()

    template = template.replace(INTEGRATION_TESTS_DIR, integration_tests_dir)
    template = template.replace(TEST_DATA_DIR, test_data_dir)
    template = template.replace(RESOURCES_DIR, resources_dir)

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

class HailTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):

        # render test config from the template
        test_data_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        test_data_path = os.path.join(test_data_path, 'control_set_small')

        resources_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        resources_path = os.path.join(resources_path, 'resources')

        render_config('inputs_test_template.yaml', test_data_path, resources_path)

        # set up path to test config
        smoke_test_dir_path = os.path.dirname(os.path.realpath(__file__))
        test_config_path = os.path.join(smoke_test_dir_path, 'inputs_test_rendered.yaml')
        os.environ['WES_CONFIG'] = test_config_path

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

if __name__ == '__main__':
    unittest.main()

