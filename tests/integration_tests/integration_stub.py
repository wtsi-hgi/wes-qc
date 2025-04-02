import argparse
import importlib
import os
import unittest
import pytest
from typing import Optional
from unittest.mock import patch

from wes_qc import teszt

# /path/to/wes_qc must be in PYTHONPATH

# List with the test files must be located in the test directory in the repo
TEST_FILES_LIST = "../test_files_list_in_bucket.txt"

# variables for test config rendering
INTEGRATION_TESTS_DIR = "{INTEGRATION_TESTS_DIR}"
TEST_DATA_DIR = "{TEST_DATA_DIR}"
RESOURCES_DIR = "{RESOURCES_DIR}"
METADATA_DIR = "{METADATA_DIR}"
TRAINING_SETS_DIR = "{TRAINING_SETS_DIR}"
VARIANT_QC_RANDOM_FOREST_DIR = "{VARIANT_QC_RANDOM_FOREST_DIR}"
PEDIGREE_FILE_NAME = "{PEDIGREE_FILE_NAME}"

# configuration file template for the integration tests
INTEGRATION_TESTS_CONFIG_TEMPLATE = "config_test_template.yaml"
INTEGRATION_TESTS_CONFIG_RENDERED_SAVEFILE = "integration_config_rendered.yaml"


# set up explicit runhash parameter for reproducible RF training
RF_RUN_TEST_HASH = "testhash"  # manually set rf run id


def render_config(
    path_to_template: str,
    test_data_dir: Optional[str],
    resources_dir: Optional[str],
    metadata_dir: Optional[str],
    training_sets_dir: Optional[str],
    variant_qc_random_forest_dir: Optional[str],
    pedigree_file_name: Optional[str],
    savefile: str = "inputs_test_rendered.yaml",
):
    """
    Read the config template and fill in the paths.
    """
    integration_tests_dir = os.path.dirname(os.path.abspath(__file__))

    # TODO agree on test and resources folder naming convention
    test_data_dir = test_data_dir if test_data_dir else os.path.join(integration_tests_dir, "control_set")
    resources_dir = resources_dir if resources_dir else os.path.join(integration_tests_dir, "resources")
    metadata_dir = metadata_dir if metadata_dir else os.path.join(integration_tests_dir, "metadata")
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
    template = template.replace(METADATA_DIR, metadata_dir)
    template = template.replace(TRAINING_SETS_DIR, training_sets_dir)
    template = template.replace(VARIANT_QC_RANDOM_FOREST_DIR, variant_qc_random_forest_dir)
    template = template.replace(PEDIGREE_FILE_NAME, pedigree_file_name)

    with open(savefile, "w") as f:
        f.write(template)


# keep an eye on test execution order
# unittest.TestLoader.sortTestMethodsUsing = None

# DEBUG: might be useful https://spark.apache.org/docs/latest/api/python/getting_started/testing_pyspark.html
# TODO: clean up hail logs


class HailTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # get test suite path in the repo
        test_suite_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

        # specify paths to the test data and resources
        test_data_download_path = os.path.join(test_suite_path, "test_data")

        test_data_path = os.path.join(test_data_download_path, "control_set_small")
        resources_path = os.path.join(test_data_download_path, "resources")
        metadata_path = os.path.join(test_data_download_path, "metadata")
        training_sets_path = os.path.join(test_data_download_path, "training_sets")
        variant_qc_random_forest_path = os.path.join(test_data_download_path, "variant_qc_random_forest")

        print(f"Downloading data from the bucket using files list {os.path.abspath(TEST_FILES_LIST)}")
        teszt.download_test_data_using_files_list(TEST_FILES_LIST, test_data_download_path)

        integration_tests_dir = os.path.dirname(os.path.realpath(__file__))
        rendered_config_savefile = os.path.join(integration_tests_dir, INTEGRATION_TESTS_CONFIG_RENDERED_SAVEFILE)

        # Get pedigree file path from class variable, defaulting to None if not set
        pedigree_file_path = getattr(cls, "pedigree_file_path", None)

        # # render test config from the template
        render_config(
            INTEGRATION_TESTS_CONFIG_TEMPLATE,
            test_data_path,
            resources_path,
            metadata_path,
            training_sets_path,
            variant_qc_random_forest_path,
            pedigree_file_name=pedigree_file_path,
            savefile=rendered_config_savefile,
        )

        # set up path to test config
        os.environ["WES_CONFIG"] = rendered_config_savefile

    @classmethod
    def tearDownClass(cls):
        # TODO: clean up logs
        pass


class IntegrationTestsStub(HailTestCase):
    @patch(
        "argparse.ArgumentParser.parse_args",
        return_value=argparse.Namespace(
            kg_to_mt=True, kg_filter_and_prune=True, kg_pc_relate=True, kg_remove_related_samples=True, all=False
        ),
    )
    def stub_0_1_import_data(self, mock_args):
        qc_step_0_1 = importlib.import_module("0-resource_preparation.1-import_1kg")
        try:
            qc_step_0_1.main()
        except Exception as e:
            self.fail(f"Step 0.1 failed with an exception: {e}")

    def stub_0_2_import_data(self):
        qc_step_0_2 = importlib.import_module("0-resource_preparation.2-generate-truthset-ht")
        try:
            qc_step_0_2.main()
        except Exception as e:
            self.fail(f"Step 0.2 failed with an exception: {e}")

    @pytest.mark.skip(
        reason="The test depends on gnomAD table, not present in downloaded test data.\n"
        "Instead we download the resulting table from the bucket and use it as a resource.\n"
    )
    def stub_0_3_import_data(self):
        qc_step_0_3 = importlib.import_module("0-resource_preparation.3-prepare-gnomAD-table")
        try:
            qc_step_0_3.main()
        except Exception as e:
            self.fail(f"Step 0.3 failed with an exception: {e}")

    def stub_1_1_import_data(self):
        qc_step_1_1 = importlib.import_module("1-import_data.1-import_gatk_vcfs_to_hail")
        try:
            qc_step_1_1.main()
        except Exception as e:
            self.fail(f"Step 1.1 failed with an exception: {e}")

    def stub_1_2_import_data(self):
        qc_step_1_2 = importlib.import_module("1-import_data.2-import_annotations")
        try:
            qc_step_1_2.main()
        except Exception as e:
            self.fail(f"Step 1.2 failed with an exception: {e}")

    def stub_1_3_import_data(self):
        qc_step_1_3 = importlib.import_module("1-import_data.3-validate-gtcheck")
        try:
            qc_step_1_3.main()
        except Exception as e:
            self.fail(f"Step 1.3 failed with an exception: {e}")

    def stub_1_4_import_data(self):
        qc_step_1_4 = importlib.import_module("1-import_data.4-mutation-spectra_preqc")
        try:
            qc_step_1_4.main()
        except Exception as e:
            self.fail(f"Step 1_4 failed with an exception: {e}")

    ### === Sample QC === ###
    def stub_2_1_sample_qc(self):
        qc_step_2_1 = importlib.import_module("2-sample_qc.1-hard_filters_sex_annotation")
        try:
            qc_step_2_1.main()
        except Exception as e:
            self.fail(f"Step 2.1 failed with an exception: {e}")

    def stub_2_2_sample_qc(self):
        qc_step_2_2 = importlib.import_module("2-sample_qc.2-prune_related_samples")

        try:
            qc_step_2_2.main()
        except Exception as e:
            self.fail(f"Step 2.2 failed with an exception: {e}")

    @patch(
        "argparse.ArgumentParser.parse_args",
        return_value=argparse.Namespace(
            merge_and_ldprune=True, pca=True, pca_plot=True, assign_pops=True, pca_plot_assigned=True, all=False
        ),
    )
    def stub_2_3_sample_qc(self, mock_args):
        qc_step_2_3 = importlib.import_module("2-sample_qc.3-population_pca_prediction")

        try:
            qc_step_2_3.main()
        except Exception as e:
            self.fail(f"Step 2.3 failed with an exception: {e}")

    def stub_2_4_sample_qc(self):
        qc_step_2_4 = importlib.import_module("2-sample_qc.4-find_population_outliers")

        try:
            qc_step_2_4.main()
        except Exception as e:
            self.fail(f"Step 2.4 failed with an exception: {e}")

    def stub_2_5_sample_qc(self):
        qc_step_2_5 = importlib.import_module("2-sample_qc.5-filter_fail_sample_qc")
        try:
            qc_step_2_5.main()
        except Exception as e:
            self.fail(f"Step 2.5 failed with an exception: {e}")

    ### === Variant QC === ###
    # mock cli arguments
    @patch(
        "argparse.ArgumentParser.parse_args",
        return_value=argparse.Namespace(all=False, split_qc=True, trios_stats=True, inbreeding=True),
    )
    def stub_3_1_variant_qc(self, mock_args):
        qc_step_3_1 = importlib.import_module("3-variant_qc.1-split_and_family_annotate")
        try:
            qc_step_3_1.main()
        except Exception as e:
            self.fail(f"Step 3.1 failed with an exception: {e}")

    def stub_3_2_variant_qc(self):
        qc_step_3_2 = importlib.import_module("3-variant_qc.2-create_rf_ht")
        try:
            qc_step_3_2.main()
        except Exception as e:
            self.fail(f"Step 3.2 failed with an exception: {e}")

    @patch("argparse.ArgumentParser.parse_args", return_value=argparse.Namespace(manual_model_id=RF_RUN_TEST_HASH))
    def stub_3_3_variant_qc(self, mock_args):
        qc_step_3_3 = importlib.import_module("3-variant_qc.3-train_rf")
        try:
            qc_step_3_3.main()
        except Exception as e:
            self.fail(f"Step 3.3 failed with an exception: {e}")

    # mock cli arguments
    def stub_3_4_variant_qc(self):
        qc_step_3_4 = importlib.import_module("3-variant_qc.4-apply_rf")
        try:
            qc_step_3_4.main()
        except Exception as e:
            self.fail(f"Step 3.4 failed with an exception: {e}")

    # mock cli arguments
    def stub_3_5_variant_qc(self):
        qc_step_3_5 = importlib.import_module("3-variant_qc.5-annotate_ht_after_rf")
        try:
            qc_step_3_5.main()
        except Exception as e:
            self.fail(f"Step 3.5 failed with an exception: {e}")

    # mock cli arguments
    def stub_3_6_variant_qc(self):
        qc_step_3_6 = importlib.import_module("3-variant_qc.6-rank_and_bin")
        try:
            qc_step_3_6.main()
        except Exception as e:
            self.fail(f"Step 3.6 failed with an exception: {e}")

    # mock cli arguments
    def stub_3_7_variant_qc(self):
        qc_step_3_7 = importlib.import_module("3-variant_qc.7-plot_rf_output")
        try:
            qc_step_3_7.main()
        except Exception as e:
            self.fail(f"Step 3.7 failed with an exception: {e}")

    # mock cli arguments
    @patch(
        "argparse.ArgumentParser.parse_args",
        return_value=argparse.Namespace(snv=92, indel=68),
    )
    def stub_3_8_variant_qc(self, mock_args):
        qc_step_3_8 = importlib.import_module("3-variant_qc.8-select_thresholds")
        try:
            qc_step_3_8.main()
        except Exception as e:
            self.fail(f"Step 3.8 failed with an exception: {e}")

    # mock cli arguments
    @patch(
        "argparse.ArgumentParser.parse_args",
        return_value=argparse.Namespace(snv=84, indel=60),
    )
    def stub_3_9_variant_qc(self, mock_args):
        qc_step_3_9 = importlib.import_module("3-variant_qc.9-filter_mt_after_variant_qc")
        try:
            qc_step_3_9.main()
        except Exception as e:
            self.fail(f"Step 3.9 failed with an exception: {e}")

    ### === Genotype QC === ###
    @patch(
        "argparse.ArgumentParser.parse_args",
        return_value=argparse.Namespace(prepare=True, evaluate_snv=True, evaluate_indel=True, plot=True, all=False),
    )
    def stub_4_1_genotype_qc(self, mock_args):
        qc_step_4_1 = importlib.import_module("4-genotype_qc.1-compare_hard_filter_combinations")
        try:
            qc_step_4_1.main()
        except Exception as e:
            self.fail(f"Step 4.1 failed with an exception: {e}")

    def stub_4_2_genotype_qc(self):
        qc_step_4_2 = importlib.import_module("4-genotype_qc.2-apply_range_of_hard_filters")
        try:
            qc_step_4_2.main()
        except Exception as e:
            self.fail(f"Step 4.2 failed with an exception: {e}")

    def stub_4_3a_genotype_qc(self):
        qc_step_4_3a = importlib.import_module("4-genotype_qc.3a-export_vcfs_range_of_hard_filters")
        try:
            qc_step_4_3a.main()
        except Exception as e:
            self.fail(f"Step 4.3a failed with an exception: {e}")

    def stub_4_3b_genotype_qc(self):
        qc_step_4_3b = importlib.import_module("4-genotype_qc.3b-export_vcfs_stringent_filters")
        try:
            qc_step_4_3b.main()
        except Exception as e:
            self.fail(f"Step 4.3b failed with an exception: {e}")

    def stub_4_4_genotype_qc(self):
        qc_step_4_4 = importlib.import_module("4-genotype_qc.4-counts_per_sample")
        try:
            qc_step_4_4.main()
        except Exception as e:
            self.fail(f"Step 4.4 failed with an exception: {e}")

    def stub_4_5_genotype_qc(self):
        qc_step_4_5 = importlib.import_module("4-genotype_qc.5-mutation-spectra_afterqc")
        try:
            qc_step_4_5.main()
        except Exception as e:
            self.fail(f"Step 4.5 failed with an exception: {e}")
