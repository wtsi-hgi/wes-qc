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
    anndir = config["general"]["annotation_dir"]

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

    mt_out_file = path_spark(outfile)
    print(f"=== Saving as hail mt to {mt_out_file}")
    mt.write(mt_out_file, overwrite=True)
    vars, samples = mt.count()
    print(f"=== Loaded dataset: {samples} samples, {vars} variants")

    print("=== Checking for duplications ===")
    duplicated_samples = find_duplicated_samples(mt)
    duplicated_samples_count = duplicated_samples.count()
    if duplicated_samples_count > 0:
        dups_path = f"{anndir}/duplicated_samples.tsv"
        print(f"=== WARNING !!! Found {duplicated_samples_count} duplicated samples. Exporting to: {dups_path}")
        duplicated_samples.export(path_spark(dups_path))
    else:
        print("=== No duplicated samples have been found.")

    duplicated_variants = find_duplicated_variants(mt)
    duplicated_variants_count = duplicated_variants.count()
    if duplicated_variants_count > 0:
        dups_path = f"{anndir}/duplicated_variants.vcf.bgz"
        print(f"=== WARNING !!! Found {duplicated_variants_count} duplicated variants. Exporting to: {dups_path}")
        hl.export_vcf(duplicated_variants, path_spark(dups_path))
    else:
        print("=== No duplicated variants have been found.")


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


def find_duplicated_samples(mt: hl.MatrixTable) -> hl.Table:
    """
    Duplicate sample IDs from a Hail MatrixTable and return as a Hail Table.

    """

    samples = mt.cols()
    grouped = samples.group_by(samples.s).aggregate(s_count=hl.agg.count())

    # Filter for duplicates
    duplicates = grouped.filter(grouped.s_count > 1)
    return duplicates


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
