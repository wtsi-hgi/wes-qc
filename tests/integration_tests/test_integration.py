import os
import argparse
import unittest
import importlib

from typing import Optional
from unittest.mock import patch
from utils.utils import download_test_data_using_files_list


# list with the test files must be located in the test directory in the repo
TEST_FILES_LIST = "../test_files_list_in_bucket.txt"

# configuration file template for the integration tests
INTEGRATION_TESTS_CONFIG_TEMPLATE = "new_config_test_template.yaml"
INTEGRATION_TESTS_CONFIG_RENDERED_SAVEFILE = "integration_config_rendered.yaml"

# set up explicit runhash parameter for reproducible RF training
RF_RUN_TEST_HASH = "testhash"  # manually set rf run id

# variables for test config rendering
INTEGRATION_TESTS_DIR = "{INTEGRATION_TESTS_DIR}"
TEST_DATA_DIR = "{TEST_DATA_DIR}"
RESOURCES_DIR = "{RESOURCES_DIR}"
TRAINING_SETS_DIR = "{TRAINING_SETS_DIR}"
VARIANT_QC_RANDOM_FOREST_DIR = "{VARIANT_QC_RANDOM_FOREST_DIR}"


def render_config(
    path_to_template: str,
    test_data_dir: Optional[str],
    resources_dir: Optional[str],
    training_sets_dir: Optional[str],
    variant_qc_random_forest_dir: Optional[str],
    savefile: str = "inputs_test_rendered.yaml",
):
    """
    Read the config template and fill in the paths.
    """
    integration_tests_dir = os.path.dirname(os.path.abspath(__file__))

    # TODO agree on test and resources folder naming convention
    test_data_dir = test_data_dir if test_data_dir else os.path.join(integration_tests_dir, "control_set")
    resources_dir = resources_dir if resources_dir else os.path.join(integration_tests_dir, "resources")
    training_sets_dir = training_sets_dir if training_sets_dir else os.path.join(integration_tests_dir, "training_sets")
    variant_qc_random_forest_dir = (
        variant_qc_random_forest_dir
        if variant_qc_random_forest_dir
        else os.path.join(integration_tests_dir, "variant_qc_random_forest_dir")
    )

    with open(path_to_template, "r") as f:
        template = f.read()

    # TODO: make versatile
    template = template.replace(INTEGRATION_TESTS_DIR, integration_tests_dir)
    template = template.replace(TEST_DATA_DIR, test_data_dir)
    template = template.replace(RESOURCES_DIR, resources_dir)
    template = template.replace(TRAINING_SETS_DIR, training_sets_dir)
    template = template.replace(VARIANT_QC_RANDOM_FOREST_DIR, variant_qc_random_forest_dir)

    with open(savefile, "w") as f:
        f.write(template)


# keep an eye on test execution order
# unittest.TestLoader.sortTestMethodsUsing = None

# DEBUG: might be useful https://spark.apache.org/docs/latest/api/python/getting_started/testing_pyspark.html
# TODO: clean up hail logs

# /path/to/wes_qc must be in PYTHONPATH
qc_step_1_1 = importlib.import_module("1-import_data.1-import_gatk_vcfs_to_hail")
qc_step_1_4 = importlib.import_module("1-import_data.4-import_1kg")

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
        # get test suite path in the repo
        test_suite_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

        # specify paths to the test data and resources
        test_data_download_path = os.path.join(test_suite_path, "test_data")

        test_data_path = os.path.join(test_data_download_path, "control_set_small")
        resources_path = os.path.join(test_data_download_path, "resources")
        training_sets_path = os.path.join(test_data_download_path, "training_sets")
        variant_qc_random_forest_path = os.path.join(test_data_download_path, "variant_qc_random_forest")

        ref_data_path = os.path.join(
            test_data_download_path, "unit_tests", "reference_output_data"
        )  # not needed for integration tests

        print(f"Downloading data from the bucket using files list {os.path.abspath(TEST_FILES_LIST)}")
        download_test_data_using_files_list(TEST_FILES_LIST, test_data_download_path)

        integration_tests_dir = os.path.dirname(os.path.realpath(__file__))
        rendered_config_savefile = os.path.join(integration_tests_dir, INTEGRATION_TESTS_CONFIG_RENDERED_SAVEFILE)

        # # render test config from the template
        # render_config('inputs_test_template.yaml', test_data_path, resources_path) # TODO: make configurable
        render_config(
            INTEGRATION_TESTS_CONFIG_TEMPLATE,
            test_data_path,
            resources_path,
            training_sets_path,
            variant_qc_random_forest_path,
            savefile=rendered_config_savefile,
        )

        # set up path to test config
        os.environ["WES_CONFIG"] = rendered_config_savefile

    @classmethod
    def tearDownClass(cls):
        # TODO: clean up logs
        pass


class IntegrationTests(HailTestCase):
    def test_1_1_import_data(self):
        try:
            qc_step_1_1.main()
        except Exception as e:
            self.fail(f"Step 1.1 failed with an exception: {e}")

    @patch(
        "argparse.ArgumentParser.parse_args",
        return_value=argparse.Namespace(
            kg_to_mt=True, kg_filter_and_prune=True, kg_pc_relate=False, kg_remove_related_samples=False
        ),
    )
    def test_1_4_import_data(self, mock_args):
        try:
            qc_step_1_4.main()
        except Exception as e:
            self.fail(f"Step 1.4 failed with an exception: {e}")

    def test_2_1_sample_qc(self):
        try:
            qc_step_2_1.main()
        except Exception as e:
            self.fail(f"Step 2.1 failed with an exception: {e}")

    def test_2_2_sample_qc(self):
        try:
            qc_step_2_2.main()
        except Exception as e:
            self.fail(f"Step 2.2 failed with an exception: {e}")

    # mock cli arguments
    @patch(
        "argparse.ArgumentParser.parse_args",
        return_value=argparse.Namespace(kg_to_mt=True, run=True, merge=True, filter=True, pca=True, assign_pops=True),
    )
    def test_2_3_sample_qc(self, mock_args):
        try:
            qc_step_2_3.main()
        except Exception as e:
            self.fail(f"Step 2.3 failed with an exception: {e}")

    def test_2_4_sample_qc(self):
        try:
            qc_step_2_4.main()
        except Exception as e:
            self.fail(f"Step 2.4 failed with an exception: {e}")

    def test_2_5_sample_qc(self):
        try:
            qc_step_2_5.main()
        except Exception as e:
            self.fail(f"Step 2.5 failed with an exception: {e}")

    # mock cli arguments
    @patch("argparse.ArgumentParser.parse_args", return_value=argparse.Namespace(all=True, truth=True, annotation=True))
    def test_3_1_variant_qc(self, mock_args):
        try:
            qc_step_3_1.main()
        except Exception as e:
            self.fail(f"Step 3.1 failed with an exception: {e}")

    def test_3_2_variant_qc(self):
        try:
            qc_step_3_2.main()
        except Exception as e:
            self.fail(f"Step 3.2 failed with an exception: {e}")

    @patch("argparse.ArgumentParser.parse_args", return_value=argparse.Namespace(manual_runhash=RF_RUN_TEST_HASH))
    def test_3_3_variant_qc(self, mock_args):
        try:
            qc_step_3_3.main()
        except Exception as e:
            self.fail(f"Step 3.3 failed with an exception: {e}")

    # mock cli arguments
    @patch("argparse.ArgumentParser.parse_args", return_value=argparse.Namespace(runhash=RF_RUN_TEST_HASH))
    def test_3_4_variant_qc(self, mock_args):
        # set up correct PYSPARK_PYTHON env variable
        try:
            qc_step_3_4.main()
        except Exception as e:
            self.fail(f"Step 3.4 failed with an exception: {e}")

    # mock cli arguments
    @patch("argparse.ArgumentParser.parse_args", return_value=argparse.Namespace(runhash=RF_RUN_TEST_HASH))
    def test_3_5_variant_qc(self, mock_args):
        try:
            qc_step_3_5.main()
        except Exception as e:
            self.fail(f"Step 3.5 failed with an exception: {e}")

    # mock cli arguments
    @patch("argparse.ArgumentParser.parse_args", return_value=argparse.Namespace(runhash=RF_RUN_TEST_HASH))
    def test_3_6_variant_qc(self, mock_args):
        try:
            qc_step_3_6.main()
        except Exception as e:
            self.fail(f"Step 3.6 failed with an exception: {e}")

    # mock cli arguments
    @patch("argparse.ArgumentParser.parse_args", return_value=argparse.Namespace(runhash=RF_RUN_TEST_HASH))
    def test_3_7_variant_qc(self, mock_args):
        try:
            qc_step_3_7.main()
        except Exception as e:
            self.fail(f"Step 3.7 failed with an exception: {e}")

    # mock cli arguments
    @patch(
        "argparse.ArgumentParser.parse_args",
        return_value=argparse.Namespace(runhash=RF_RUN_TEST_HASH, snv=92, indel=68),
    )
    def test_3_8_variant_qc(self, mock_args):
        try:
            qc_step_3_8.main()
        except Exception as e:
            self.fail(f"Step 3.8 failed with an exception: {e}")

    # mock cli arguments
    @patch(
        "argparse.ArgumentParser.parse_args",
        return_value=argparse.Namespace(runhash=RF_RUN_TEST_HASH, snv=84, indel=60),
    )
    def test_3_9_variant_qc(self, mock_args):
        try:
            qc_step_3_9.main()
        except Exception as e:
            self.fail(f"Step 3.9 failed with an exception: {e}")

    # mock cli arguments
    @patch("argparse.ArgumentParser.parse_args", return_value=argparse.Namespace(dp=5, gq=10, ab=0.2))
    def test_4_1_genotype_qc(self, mock_args):
        try:
            qc_step_4_1.main()
        except Exception as e:
            self.fail(f"Step 4.1 failed with an exception: {e}")

    # mock cli arguments
    @patch("argparse.ArgumentParser.parse_args", return_value=argparse.Namespace(runhash=RF_RUN_TEST_HASH))
    def test_4_1a_genotype_qc(self, mock_args):
        try:
            qc_step_4_1a.main()
        except Exception as e:
            self.fail(f"Step 4.1a failed with an exception: {e}")

    def test_4_2_genotype_qc(self):
        try:
            qc_step_4_2.main()
        except Exception as e:
            self.fail(f"Step 4.2 failed with an exception: {e}")

    def test_4_3_genotype_qc(self):
        try:
            qc_step_4_3.main()
        except Exception as e:
            self.fail(f"Step 4.3 failed with an exception: {e}")

    # mock cli arguments
    @patch("argparse.ArgumentParser.parse_args", return_value=argparse.Namespace(runhash=RF_RUN_TEST_HASH))
    def test_4_3a_genotype_qc(self, mock_args):
        try:
            qc_step_4_3a.main()
        except Exception as e:
            self.fail(f"Step 4.3a failed with an exception: {e}")

    # mock cli arguments
    @patch("argparse.ArgumentParser.parse_args", return_value=argparse.Namespace(runhash=RF_RUN_TEST_HASH))
    def test_4_3b_genotype_qc(self, mock_args):
        try:
            qc_step_4_3b.main()
        except Exception as e:
            self.fail(f"Step 4.3b failed with an exception: {e}")


if __name__ == "__main__":
    unittest.main()
