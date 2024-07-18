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

        # set up path to test config
        smoke_test_dir_path = os.path.dirname(os.path.realpath(__file__))
        test_config_path = os.path.join(smoke_test_dir_path, 'inputs_test.yaml')
        os.environ['WES_CONFIG'] = test_config_path

    @classmethod
    def tearDownClass(cls):
        pass


class TestImportData(HailTestCase):
    def test_all(self):
        try:
            qc_step_1_1.main()
            hl.stop()
            qc_step_2_1.main()
            hl.stop()
            qc_step_2_2.main()
            hl.stop()
        except Exception as e:
            self.fail(f'QC failed with an exception: {e}')

        # no errors
        # self.assertEqual((self.output_mt_path, hfs.is_dir(self.output_mt_path)), 
        #                  (self.output_mt_path, True))

if __name__ == '__main__':
    unittest.main()

# if __name__ == '__main__':
#     unittest.main(argv=[''], verbosity=0, exit=False)

