import os
import unittest
import importlib
import hail as hl
import hailtop.fs as hfs
from pyspark import SparkContext

# DEBUG: might be useful https://spark.apache.org/docs/latest/api/python/getting_started/testing_pyspark.html

# TODO: don't generate logs, or clean them up in the tearDown() method (e.g. using glob())

# TODO: find on the internet how to run the entire test suite using a custom cli command 
# (which is spark-submit in our case)

# /path/to/wes_qc must be in PYTHONPATH
qc_step_1 = importlib.import_module("1-import_data.1-import_gatk_vcfs_to_hail")

class TestImportData(unittest.TestCase):
    def setUp(self):
        # TODO: try creating simulated test vcf data here

        # the data is located one dir up
        test_dataset_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__))) # TODO: move to class variables?
        self.import_vcf_dir = f"file://{os.path.join(test_dataset_path, 'control_set_small')}"

        # the output dir will be in the current dir. TODO: think of better design
        test_outdir_path = os.path.dirname(os.path.realpath(__file__))
        self.mtdir = f"file://{os.path.join(test_outdir_path, 'matrixtables_test')}"
        self.tmp_dir = f"file://{os.path.join(test_outdir_path, 'tmp_test')}"
        self.vcf_header = '' # not present for test dataset. TODO: create one to ensure test completeness

    def test_import(self):
        # TODO: create spark session during class initialisation
        sc = SparkContext.getOrCreate()
        hadoop_config = sc._jsc.hadoopConfiguration()
        hl.init(sc=sc, tmp_dir=self.tmp_dir, default_reference="GRCh38", idempotent=True)

        # load VCFs
        qc_step_1.load_vcfs_to_mt(self.import_vcf_dir, self.mtdir, self.tmp_dir, self.vcf_header)

        # the name of the output .mt file is  "gatk_unprocessed.mt"
        # hardcoded in 1-import_gatk_vcfs_to_hail.py:load_vcfs_to_mt()
        # TODO: make output filename variable and refactor test
        output_mt_filename = 'gatk_unprocessed.mt'
        output_mt_path = os.path.join(self.mtdir, output_mt_filename)

        # you can do smth like that to test if the directory was created
        self.assertEqual((output_mt_path, hfs.is_dir(output_mt_path)), (output_mt_path, True))

    def tearDown(self):
        # remove temporary test data
        hfs.rmtree(self.mtdir)
        hfs.rmtree(self.tmp_dir)

if __name__ == '__main__':
    unittest.main(argv=[''], verbosity=0, exit=False)
