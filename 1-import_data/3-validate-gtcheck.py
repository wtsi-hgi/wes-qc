"""
This module provides functionality to validate and analyze the consistency between
whole-exome sequencing (WES) and microarray genotyping data.

**Definitions**
* WES - whole-exome or whole-genome sequencing
* Array - data from array-based genotyping
* Mapping file - the table contains pairs between WES ID and MicroARRAY ID

The issue: Very often, the mapping file quality is bad.
It contains non-unique IDs (for example, one WES IDs pairing with multiple ARRAY IDs,
and multiple mismatches (for example, IDs present in the real data but not covered by mapping)

**Definition of ID sources**
* Data IDs - IDs, present in the real data: VCFs for WES, VCFs or FAM for ARRAY
* Mapping IDs -- IDs, exiting in the mapping file

**Required input data**
* List of WES IDs from the data - extracted from the input Hail matrixtable
* List of ARRAY IDs from data
* Mapping file

The mapping file should contain at least two columns - one for WES ID, and one for ARRAY ID
If you have mapping information in any other form, you need to manually convert it.


"""

import os

import pandas as pd
from collections import Counter
from pathlib import Path
import hail as hl

from utils.config import path_spark
from utils.utils import parse_config
from wes_qc import hail_utils

Id = str
IdList = list[Id]
IdSet = set[Id]

# ======== The first part - validation of the WES-microarray mapping file =======


def check_id_uniquness(ids: IdList) -> tuple[IdSet, IdSet]:
    """
    Returns unique and non-unique IDs from the list of IDs
    """

    counts = Counter(ids)
    unique = set((key for key, count in counts.items() if count == 1))
    non_unique = set((key for key, count in counts.items() if count > 1))
    return unique, non_unique


def dump_ids(ids_list: IdList, filename: Path):
    """
    Dumps to the text file all items from the list
    """
    with open(filename, "w") as f:
        f.write("\n".join(ids_list))


def report_id_uniquness(unique, non_unique, data_type, dump_file_name=None):
    print(f"{len(unique)} unique IDs in {data_type}")
    if len(non_unique) > 0:
        print(f"{len(non_unique)} non-unique IDs in {data_type}")
        if dump_file_name is not None:
            print(f"Non-unique samples dumped to {dump_file_name}")
            dump_ids(non_unique, dump_file_name)
    else:
        print(f"All IDs in {data_type} are unique")


def check_ids_consistency(data_ids, map_ids, caption: str, dump_prefix: Path):
    data_ids_unique, data_ids_non_unique = check_id_uniquness(data_ids)
    map_ids_unique, map_ids_non_unique = check_id_uniquness(map_ids)

    data_ids, map_ids = set(data_ids), set(map_ids)

    data_no_map = data_ids - map_ids
    if len(data_no_map) > 0:
        data_no_map_name = dump_prefix + ".ids_data_not_in_map.txt"
        print(f"{len(data_no_map)} {caption} IDs from data not present in mapping: {data_no_map_name}")
        print(
            f"   --> Including {len(data_no_map & data_ids_unique)} unique IDs, and {len(data_no_map & data_ids_non_unique)} non-uniqie IDS"
        )
        dump_ids(data_no_map, data_no_map_name)

    else:
        print(f"All {caption} data IDs present in mapping table")
    map_no_data = map_ids - data_ids
    if len(map_no_data) > 0:
        map_no_data_name = dump_prefix + ".ids_map_not_in_data.txt"
        print(f"{len(map_no_data)} {caption} IDs from map not present in data: {map_no_data_name}")
        print(
            f"   --> Including {len(map_no_data & map_ids_unique)} unique IDs, and {len(map_no_data & map_ids_non_unique)} non-uniqie IDS"
        )
        dump_ids(map_no_data, map_no_data_name)
    else:
        print(f"All {caption} IDs from mapping exist in data")
    return data_no_map, map_no_data


def mapping_mark_duplicates(mapping: pd.DataFrame, wes_id_col: str = "wes_id", microarray_id_col: str = "array_id"):
    mapping = mapping.copy()
    mapping["duplicated_wes"] = mapping[wes_id_col].duplicated(keep=False)
    mapping["duplicated_array"] = mapping[microarray_id_col].duplicated(keep=False)
    mapping["duplicated_id_any"] = mapping.duplicated_wes | mapping.duplicated_array
    mapping["duplicated_id_both"] = mapping.duplicated_wes & mapping.duplicated_array
    return mapping


def report_id_stats(ids: IdList, caption: str, stats_non_unique_file: Path) -> None:
    print(f"\n=== {caption} samples validation ===")
    print(f"Total {len(ids)} IDs in {caption} data")
    ids_unique, ids_non_unique = check_id_uniquness(ids)
    report_id_uniquness(ids_unique, ids_non_unique, caption, stats_non_unique_file)


def report_mapping_duplicated_stats(mapping, duplicated_samples_name) -> None:
    print("\n=== Marking duplicates in ID mapping ===")
    if mapping.duplicated_id_any.any():
        mapping_dubl = mapping[mapping.duplicated_id_any]
        print(f"Identified {len(mapping_dubl)} samples with non-unique WES<->array mapping: {duplicated_samples_name}")
        mapping_dubl.to_csv(duplicated_samples_name)
    else:
        print("WES<->array mapping is unique")


def verify_wes2array_mapping(
    ids_wes_data,
    ids_microarray_data,
    mapping,
    wes_id_col="wes_id",
    microarray_id_col="microarray_id",
    results_file_prefix="",
):
    print("\n=== Mapping stats ===")
    print(f"Total: {len(mapping)} records in mapping")

    ids_wes_mapping = mapping[wes_id_col].to_list()
    ids_microarray_mapping = mapping[microarray_id_col].to_list()

    ids_wes_mapping_unique, ids_wes_mapping_non_unique = check_id_uniquness(ids_wes_mapping)
    report_id_uniquness(
        ids_wes_mapping_unique,
        ids_wes_mapping_non_unique,
        "WES mapping",
        results_file_prefix + ".WES_mapping_non-unique.txt",
    )

    ids_microarray_mapping_unique, ids_microarray_mapping_non_unique = check_id_uniquness(ids_microarray_mapping)
    report_id_uniquness(
        ids_microarray_mapping_unique,
        ids_microarray_mapping_non_unique,
        "ARAY mapping",
        results_file_prefix + ".ARRAY_mapping_non-unique.txt",
    )

    print("\n=== Mapping consistency checking ===")
    check_ids_consistency(ids_wes_data, ids_wes_mapping, "WES", results_file_prefix + ".WES")
    check_ids_consistency(ids_microarray_data, ids_microarray_mapping, "ARRAY", results_file_prefix + ".ARRAY")

    mapping_marked_duplicated = mapping_mark_duplicates(mapping, wes_id_col, microarray_id_col)
    report_mapping_duplicated_stats(mapping_marked_duplicated, results_file_prefix + ".non-unique-pairs.csv")


def mapping_clean_missing(
    wes2array_mapping: pd.DataFrame,
    ids_microarray_data: IdList,
    microarray_id_col="microarray_id",
) -> pd.DataFrame:
    print("Removing from the mapping all microarray IDs not present in the real microarray data")
    ids_to_keep = wes2array_mapping[microarray_id_col].isin(ids_microarray_data)
    wes2array_clean_missing = wes2array_mapping.loc[ids_to_keep, :]
    ids_array_mapping_corrected_clean_missing = set(wes2array_clean_missing[microarray_id_col])
    print("Survived unique Array IDs: ", len(ids_array_mapping_corrected_clean_missing))
    return wes2array_clean_missing


# ======== The second part - validation of the gtcheck data consistency =======


def main() -> None:
    # = STEP SETUP = #
    config = parse_config()

    tmp_dir = config["general"]["tmp_dir"]

    # = STEP PARAMETERS = #
    dataset = config["general"]["dataset_name"]
    wes_id_col = config["step1"]["validate_gtcheck"]["wes_id_col"]
    microarray_id_col = config["step1"]["validate_gtcheck"]["microarray_id_col"]

    # = STEP DEPENDENCIES = #
    mtpath = config["step1"]["validate_verifybamid"]["mt_metadata_annotated"]
    wes_microarray_mapping = config["step1"]["validate_gtcheck"]["wes_microarray_mapping"]
    microarray_ids = config["step1"]["validate_gtcheck"]["microarray_ids"]

    # = STEP OUTPUTS = #
    gtcheck_results_folder = config["step1"]["validate_gtcheck"]["gtcheck_results_folder"]
    gtcheck_results_prefix = os.path.join(gtcheck_results_folder, dataset)

    # = STEP LOGIC = #
    if wes_microarray_mapping is None:
        exit()  # TODO: mock_output matrixtable

    # === The first part - validation if the  WES-micro array mapping consistency ===

    Path(gtcheck_results_folder).mkdir(parents=True, exist_ok=True)

    # TODO: changed to local data loading for debug purposes. Uncomment before running on cluster.

    _ = hail_utils.init_hl(tmp_dir)

    # Loading the list of WES samples from Hial matrixtable
    mt = hl.read_matrix_table(path_spark(mtpath))
    ids_wes_data = mt.s.collect()
    """
    with open(
        "/lustre/scratch126/teams/hgi/gz3/wes_qc_pycharm/tests/test_data/metadata/control_set_small.wes_samples.txt"
    ) as warray:
        ids_wes_data = [s.strip() for s in warray.readlines()]
    """
    # Loading the list of microarray samples form an external file
    with open(microarray_ids) as farray:
        ids_microarray_data = [s.strip() for s in farray.readlines()]

    report_id_stats(ids_wes_data, "WES data", f"{gtcheck_results_folder}/{dataset}.non-unique-ids.data.WES.txt")
    report_id_stats(
        ids_microarray_data, "ARRAY data", f"{gtcheck_results_folder}/{dataset}.non-unique-ids.data.microarray.txt"
    )

    # - Validating mapping table
    wes2array = pd.read_csv(wes_microarray_mapping, sep="\t")
    verify_wes2array_mapping(
        ids_wes_data,
        ids_microarray_data,
        wes2array,
        wes_id_col,
        microarray_id_col,
        results_file_prefix=gtcheck_results_prefix,
    )

    # - Correcting the mapping file - Removing microarray IDs not present in the real microarray data
    wes2array = mapping_clean_missing(wes2array, ids_microarray_data, microarray_id_col)


if __name__ == "__main__":
    main()
