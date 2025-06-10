# Load GATK VCFs into hail and save as matrixtable
import re
import hail as hl
from utils.utils import parse_config, path_spark
from wes_qc import hail_utils

# DEBUG: for some reason, paths prefix is `file:`, not a `file://`
VCF_PATTERN = re.compile("file:.*vcf.b?gz$")


def load_vcfs_to_mt(config):
    """
    load VCFs and save as hail mt.
    Save mt as outdir/gatk_unprocessed.mt

    ### Config fields
    ```
    step1.gatk_vcf_header_infile
    step1.gatk_vcf_indir
    step1.gatk_mt_outfile
    ```
    """
    indir, header, outfile = (
        config["step1"]["gatk_vcf_indir"],
        config["step1"].get("gatk_vcf_header_infile"),  # optional
        config["step1"]["gatk_mt_outfile"],
    )

    objects = hl.utils.hadoop_ls(path_spark(indir))

    # get paths of all vcf files
    vcfs = [vcf["path"] for vcf in objects if VCF_PATTERN.match(vcf["path"])]
    print(f"info: Found {len(vcfs)} VCFs in {indir}")
    # create and save MT
    if header:
        print("info: Loading VCFs with header")
        mt = hl.import_vcf(vcfs, array_elements_required=False, force_bgz=True, header_file=header)
    else:
        print("info: Loading VCFs WITHOUT header")
        mt = hl.import_vcf(vcfs, array_elements_required=False, force_bgz=True)


    print("=== Checking for duplications ===")
    duplicated_variants = find_duplicated_variants(mt)
    if duplicated_variants.count() > 0:
        dups_path = "{outfile}.duplicated_variants.vcf"
        print(f"=== WARNING !!! Found duplicated variants. Exporting to: {dups_path}")
        hl.export_vcf(duplicated_variants, path_spark(dups_path))
    else:
        print("=== Check for duplicates variants: No duplicates has been found. Everything looks good...")

    mt_out_file = path_spark(outfile)
    print(f"Saving as hail mt to {mt_out_file}")
    mt.write(mt_out_file, overwrite=True)

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


def main():
    # = STEP SETUP = #
    config = parse_config()
    tmp_dir = config["general"]["tmp_dir"]

    # = STEP LOGIC = #
    _ = hail_utils.init_hl(tmp_dir)

    # load VCFs
    load_vcfs_to_mt(config)

    # sc.stop() # DEBUG


if __name__ == "__main__":
    main()
