import unittest
import importlib
import pathlib as pl
import shutil as sh

# /path/to/parent/of/wes_qc must be in PYTHONPATH
qc_step_1 = importlib.import_module("wes_qc.1-import_data.1-import_gatk_vcfs_to_hail")

class TestImportData(unittest.TestCase):
    def setUp(self):
        # create temporary test vcf data here
        self.indir = 'input_data'
        self.outdir = pl.Path('output_data/')
        self.outdir.mkdir()
        self.tmp_dir = '' # not used in load_vsf_to_mt
        self.header = '' # not used in code

    def test_import(self):

        # qc_step_1.load_vcfs_to_mt(indir, outdir, tmp_dir, header)

        # TODO: add test for matrixtable creation
        # the name of the output .mt file is  "gatk_unprocessed.mt"
        # hardcoded in 1-import_gatk_vcfs_to_hail.py:load_vcfs_to_mt()

        # TODO: make output filename variable
        output_mt_filename = 'gatk_unprocessed.mt'
        output_mt_path = self.outdir / output_mt_filename

        # you can do smth like that to test if the directory was created
        self.assertEqual((str(output_mt_path), output_mt_path.is_dir()), (str(output_mt_path), True))

    def tearDown(self):
        # remove temporary test data
        sh.rmtree(self.outdir)

if __name__ == '__main__':
    unittest.main()
