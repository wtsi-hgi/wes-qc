# Load GATK VCFs into hail and save as matrixtable
import hail as hl
import os

from utils.utils import parse_config
from wes_qc import hail_utils


def load_vcfs_to_mt(indir: str, outdir: str, header: str) -> None:
    """
    load VCFs and save as hail mt
    """
    objects = hl.utils.hadoop_ls(indir)
    print("Collecting VCFs")
    vcfs = [vcf["path"] for vcf in objects if (vcf["path"].startswith("file") and vcf["path"].endswith("vcf.gz"))]
    print("Loading VCFs")
    # create and save MT
    mt = hl.import_vcf(vcfs, array_elements_required=False, force_bgz=True, header_file=header)
    print("Saving as hail mt")
    mt_out_file = os.path.join(outdir, "gatk_unprocessed.mt")
    mt.write(mt_out_file, overwrite=True)


def main() -> None:
    # set up input variables
    inputs = parse_config()

    data_root: str = inputs["data_root"]
    dataset_name: str = inputs["dataset_name"]

    import_vcf_dir = os.path.join(data_root, inputs["gatk_import_lustre_dir"])
    mtdir: str = os.path.join(data_root, inputs["matrixtables_lustre_dir"])
    annot_dir: str = os.path.join(data_root, inputs["annotation_lustre_dir"])

    tmp_dir: str = inputs["tmp_dir"]

    vcf_header = os.path.join(data_root, annot_dir, inputs["gatk_vcf_header"])

    # initialise hail
    sc = hail_utils.init_hl(tmp_dir)

    # load VCFs
    load_vcfs_to_mt("file://" + import_vcf_dir, "file://" + mtdir, "file://" + vcf_header)
    hail_utils.stop_hl(sc)


if __name__ == "__main__":
    main()
