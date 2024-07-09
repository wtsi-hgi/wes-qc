# Load GATK VCFs into hail and save as matrixtable
import hail as hl
import os

from utils.utils import parse_config
from wes_qc import hail_utils


def load_vcfs_to_mt(indir: str, header: str) -> hl.MatrixTable:
    """
    load VCFs and save as hail mt
    """
    objects = hl.utils.hadoop_ls(indir)
    print("Collecting VCFs")
    vcfs = [vcf["path"] for vcf in objects if (vcf["path"].startswith("file") and vcf["path"].endswith("vcf.gz"))]
    print(f"Collected VCFs: {len(vcfs)}")
    print("=== List of VCFs to load:  ")
    print("\n".join(vcfs))

    # create and save MT
    mt = hl.import_vcf(vcfs, array_elements_required=False, force_bgz=True, header_file=header)
    return mt


def find_duplicated_variants(mt: hl.MatrixTable) -> hl.Table:
    # Group by chromosome, position, and reference allele
    grouped = mt.rows()
    grouped = grouped.group_by(variation=(grouped.locus, grouped.alleles[0])).aggregate(var_count=hl.agg.count())

    # Filter for duplicates
    duplicates = grouped.filter(grouped.var_count > 1)
    duplicates = duplicates.annotate(dup_locus=duplicates.variation[0])
    duplicates = duplicates.key_by("dup_locus")

    mt_duplicated_records = mt.filter_rows(hl.is_defined(duplicates[mt.locus]))
    return mt_duplicated_records.rows()


def exclude_variations(mt: hl.MatrixTable, variants_file: str) -> hl.MatrixTable:
    """
    Exclude variations from a Hail MatrixTable based on a TSV file of variations to exclude.

    Parameters:
    mt (hl.MatrixTable): Input Hail MatrixTable
    variants_file (str): Path to the VCF file containing variations to exclude.

    Returns:
    hl.MatrixTable: MatrixTable with specified variations excluded
    """
    # Read the variants from the TSV file into a Hail Table
    ht_variants = hl.import_vcf(variants_file).rows()

    # Define a filter expression for the variations to exclude
    def filter_expr(row: hl.MatrixTable) -> hl.expr:
        return hl.case().when(hl.is_defined(ht_variants.index(row.locus, row.alleles)), True).default(False)

    # Filter the MatrixTable to exclude the specified variations
    mt_filtered = mt.filter_rows(filter_expr(mt.row), keep=False)

    return mt_filtered


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
    mt = load_vcfs_to_mt("file://" + import_vcf_dir, "file://" + vcf_header)

    # Exclude variants from blacklist
    if inputs["exclude_variations_file"] != "":
        # Manual remove variations from the list
        variants_file = os.path.join(annot_dir, inputs["exclude_variations_file"])
        print(f"Excluding variants from VCF: {variants_file}")
        print(f"Initial variations: {mt.count_rows()}")
        mt = exclude_variations(mt, "file://" + variants_file)
        print(f"Variations after filtering: {mt.count_rows()}")

    print("Saving as hail mt")
    mt_out_file = "file://" + os.path.join(mtdir, "gatk_unprocessed.mt")
    mt.write(mt_out_file, overwrite=True)

    # Test for duplicates
    mt = hl.read_matrix_table("file://" + os.path.join(mtdir, "gatk_unprocessed.mt"))
    duplicated_variants = find_duplicated_variants(mt)

    if duplicated_variants.count() > 0:
        dups_path = os.path.join(annot_dir, "duplicated_variants.vcf")
        print(f"=== WARNING !!! Found duplicated variants. Exporting to {dups_path}")
        hl.export_vcf(duplicated_variants, "file://" + dups_path)
    else:
        print("=== Check for duplicates variants: No duplicates has been found. Everything looks good...")

    hail_utils.stop_hl(sc)


if __name__ == "__main__":
    main()
