import os
import unittest
import importlib
import hail as hl
import hailtop.fs as hfs
from pyspark import SparkContext

# /path/to/wes_qc/ must be in PYTHONPATH

# step 1 is needed to create mt from vcf data
# TODO: think of a better test setup

qc_step_1 = importlib.import_module("1-import_data.1-import_gatk_vcfs_to_hail")
qc_step_2 = importlib.import_module("2-sample_qc.1-hard_filters_sex_annotation")

class TestHardFiltersSexAnnotationData(unittest.TestCase):
    def setUp(self):

        # the test data is located one dir up
        test_dataset_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__))) # TODO: move to class variables?
        self.import_vcf_dir = f"file://{os.path.join(test_dataset_path, 'control_set_small')}"

        # the output dir will be in the current dir
        test_outdir_path = os.path.dirname(os.path.realpath(__file__))
        self.mtdir = f"file://{os.path.join(test_outdir_path, 'matrixtables_test')}"
        self.tmp_dir = f"file://{os.path.join(test_outdir_path, 'tmp_test')}"

        # DEBUG: in case the test function don't create the dirs itself
        # maybe will not be needed after having all tests implemented
        hfs.mkdir(self.mtdir)
        hfs.mkdir(self.tmp_dir)

        self.vcf_header = '' # not present for test dataset. TODO: create one to ensure test completeness

        # Set up spark session
        # sc = pyspark.SparkContext()
        sc = SparkContext.getOrCreate()
        hadoop_config = sc._jsc.hadoopConfiguration()
        hl.init(sc=sc, tmp_dir=self.tmp_dir, default_reference="GRCh38", idempotent=True)

        # Create mt table for tests
        qc_step_1.load_vcfs_to_mt(self.import_vcf_dir, self.mtdir, self.tmp_dir, self.vcf_header)

        output_mt_filename = 'gatk_unprocessed.mt'
        self.mt_table = os.path.join(self.mtdir, output_mt_filename)

    def test_apply_hard_filters(self):
        # read mt and apply hard filters
        mt_unfiltered = hl.read_matrix_table(self.mt_table)
        mt_filtered = qc_step_2.apply_hard_filters(mt_unfiltered, self.mtdir)

        # output filename is hardcoded in 2-sample_qc/1-hard_filters_sex_annotation.py:apply_hard_filters()
        # TODO: make output filename variable and refactor test
        filtered_mt_filename = 'mt_hard_filters_annotated.mt'
        filtered_mt_path = os.path.join(self.mtdir, filtered_mt_filename)

        self.assertEqual((filtered_mt_path, hfs.is_dir(filtered_mt_path)), (filtered_mt_path, True))

    def test_impute_sex(self):
        # TODO: implement similarly to test_apply_hard_filters
        pass

    def test_identify_inconsitencies(self):
        # TODO: implement similarly to test_apply_hard_filters
        pass

    def tearDown(self):
        # remove temporary test data
        hfs.rmtree(self.mtdir)
        hfs.rmtree(self.tmp_dir)

if __name__ == '__main__':
    unittest.main(argv=[''], verbosity=0, exit=False)
