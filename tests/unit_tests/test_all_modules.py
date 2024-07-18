import os
import unittest
import importlib
import hail as hl
import hailtop.fs as hfs
from pyspark import SparkContext

# ensure sequential test execution order
# unittest.TestLoader.sortTestMethodsUsing = None

# DEBUG: might be useful https://spark.apache.org/docs/latest/api/python/getting_started/testing_pyspark.html
# TODO: clean up hail logs

# /path/to/wes_qc must be in PYTHONPATH
qc_step_1_1 = importlib.import_module("1-import_data.1-import_gatk_vcfs_to_hail")
qc_step_2_1 = importlib.import_module("2-sample_qc.1-hard_filters_sex_annotation")
qc_step_2_2 = importlib.import_module("2-sample_qc.2-prune_related_samples")

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


class TestImportData(HailTestCase):
    @classmethod
    def setUpClass(cls):
        # set up Hail from the base class
        super(TestImportData, cls).setUpClass()

        # define parameters needed for functions to be tested

        # path to dir with the test data in the repo
        cls.test_dataset_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

        # QC Step 1.1 inputs
        cls.import_vcf_dir = f"file://{os.path.join(cls.test_dataset_path, 'control_set_small')}"
        cls.mtdir = f"file://{os.path.join(cls.test_outdir_path, 'matrixtables_test')}"
        cls.vcf_header = '' # not available for test dataset. TODO: create one to ensure test completeness
        # QC Step 1.1 outputs
        # TODO: make output filename variable and refactor the test
        cls.output_mt_path = os.path.join(cls.mtdir, 'gatk_unprocessed.mt') 

        # QC Step 2.1 outputs
        # TODO: make output filename variable and refactor the test
        cls.output_filtered_mt_path = os.path.join(cls.mtdir, 'mt_hard_filters_annotated.mt') 
        
    # QC Step 1.1
    def test_1_load_vcfs_to_mt(self):
        qc_step_1_1.load_vcfs_to_mt(self.import_vcf_dir, self.mtdir, self.tmp_dir, self.vcf_header)

        # TODO: change to comparison with reference output data
        # check MatrixTable creation
        self.assertEqual((self.output_mt_path, hfs.is_dir(self.output_mt_path)), 
                         (self.output_mt_path, True))


    # QC Step 2.1
    def test_2_apply_hard_filters(self):
        # read output of the step 1.1
        mt_unfiltered = hl.read_matrix_table(self.output_mt_path)
        mt_filtered = qc_step_2_1.apply_hard_filters(mt_unfiltered, self.mtdir)

        # TODO: change to comparison with reference output data
        self.assertEqual((self.output_filtered_mt_path, hfs.is_dir(self.output_filtered_mt_path)),
                         (self.output_filtered_mt_path, True))

    @classmethod
    def tearDownClass(cls):
        # clean up the directories created by the QC steps
        hfs.rmtree(cls.mtdir)

        super(TestImportData, cls).tearDownClass()


if __name__ == '__main__':
    unittest.main()

# if __name__ == '__main__':
#     unittest.main(argv=[''], verbosity=0, exit=False)

