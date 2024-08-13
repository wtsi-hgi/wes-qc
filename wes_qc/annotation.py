"""
A set of functions to manipulate with sample annotaions:
change sample names, sample fields, etc.
"""

from typing import Dict
import hail as hl
from wes_qc.hail_utils import trace_analysis_step


def rename_samples(mt: hl.MatrixTable, mapping_file: str) -> hl.MatrixTable:
    """
    Changes sample names according to the provided mapping_file
    """
    mt = trace_analysis_step(mt, f"Renaming table: {mapping_file}")

    # Read the mapping file into a Hail Table
    ht = hl.import_table(mapping_file, no_header=True, comment="#", delimiter="\t")

    # Rename the columns for easier access
    ht = ht.rename({"f0": "old_id", "f1": "new_id"})

    # Convert the Hail Table to a dictionary
    rename_dict: Dict[str, str] = dict(hl.tuple([ht.old_id, ht.new_id]).collect())

    # Annotate the samples with the new IDs
    mt = mt.key_cols_by(s=hl.literal(rename_dict).get(mt.s))

    return mt
