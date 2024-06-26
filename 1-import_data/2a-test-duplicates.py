# Load GATK VCFs into hail and save as matrixtable
import hail as hl
import os

from utils.utils import parse_config
from wes_qc import hail_utils
from typing import List, Tuple


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


def main() -> None:
    # set up input variables
    inputs = parse_config()

    data_root: str = inputs["data_root"]
    mtdir: str = os.path.join(data_root, inputs["matrixtables_lustre_dir"])
    annot_dir: str = os.path.join(data_root, inputs["annotation_lustre_dir"])

    tmp_dir: str = inputs["tmp_dir"]

    # initialise hail
    sc = hail_utils.init_hl(tmp_dir)

    mt = hl.read_matrix_table("file://" + os.path.join(mtdir, "gatk_unprocessed.mt"))
    duplicated_variants = find_duplicated_variants(mt)
    # duplicated_variants.export(os.path.join(annot_dir, 'duplicated_variants_list.tsv'))

    if duplicated_variants.count() > 0:
        dups_path = os.path.join(annot_dir, "duplicated_variants_list.tsv")
        print(f"=== Found duplicated variants. Exporting to {dups_path}")

        duplicated_variants.export("file://" + dups_path)
    else:
        print("No duplicated variants - all looks good")
    hail_utils.stop_hl(sc)


if __name__ == "__main__":
    main()
