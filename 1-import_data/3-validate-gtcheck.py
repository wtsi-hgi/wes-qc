"""
This module provides functionality to validate and analyze the consistency between
whole-exome sequencing (WES) and microarray genotyping data.

**Definitions**
* WES - whole-exome or whole-genome sequencing
* Microarray - data from microarray-based genotyping
* Mapping file - the table contains pairs between WES ID and MicroARRAY ID

The issue: Very often, the mapping file quality is bad.
It contains non-unique IDs (for example, one WES IDs pairing with multiple MicroARRAY IDs,
and multiple mismatches (for example, IDs present in the real data but not covered by mapping)

**Definition of ID sources**
* Data IDs - IDs, present in the real data: VCFs for WES, VCFs or FAM for MicroARRAY
* Mapping IDs -- IDs, exiting in the mapping file

**Required input data**
* List of WES IDs from the data - extracted from the input Hail matrixtable
* List of MicroARRAY IDs from data
* Mapping file

The mapping file should contain at least two columns - one for WES ID, and one for MicroARRAY ID
If you have mapping information in any other form, you need to manually convert it.


"""

import os
from typing import Optional

import pandas as pd
from collections import Counter
from pathlib import Path
import hail as hl
import numpy as np
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


def report_id_uniquness(unique: IdSet, non_unique: IdSet, data_type: str, dump_file_name: Optional[str] = None) -> None:
    print(f"{len(unique)} unique IDs in {data_type}")
    if len(non_unique) > 0:
        print(f"{len(non_unique)} non-unique IDs in {data_type}")
        if dump_file_name is not None:
            print(f"Non-unique samples dumped to {dump_file_name}")
            dump_ids(non_unique, dump_file_name)
    else:
        print(f"All IDs in {data_type} are unique")


def check_ids_consistency(data_ids: IdList, map_ids: IdList, caption: str, dump_prefix: str) -> tuple[IdSet, IdSet]:
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


def mapping_mark_duplicates(
    mapping: pd.DataFrame, wes_id_col: str = "wes_id", microarray_id_col: str = "microarray_id"
):
    mapping = mapping.copy()
    mapping["duplicated_wes"] = mapping[wes_id_col].duplicated(keep=False)
    mapping["duplicated_microarray"] = mapping[microarray_id_col].duplicated(keep=False)
    mapping["duplicated_id_any"] = mapping.duplicated_wes | mapping.duplicated_microarray
    mapping["duplicated_id_both"] = mapping.duplicated_wes & mapping.duplicated_microarray
    return mapping


def report_id_stats(ids: IdList, caption: str, stats_non_unique_file: Path) -> None:
    print(f"\n=== {caption} samples validation ===")
    print(f"Total {len(ids)} IDs in {caption} data")
    ids_unique, ids_non_unique = check_id_uniquness(ids)
    report_id_uniquness(ids_unique, ids_non_unique, caption, stats_non_unique_file)


def report_mapping_duplicated_stats(mapping: pd.DataFrame, duplicated_samples_name: str) -> None:
    print("\n=== Marking duplicates in ID mapping ===")
    if mapping.duplicated_id_any.any():
        mapping_dubl = mapping[mapping.duplicated_id_any]
        print(
            f"Identified {len(mapping_dubl)} samples with non-unique WES<->microarray mapping: {duplicated_samples_name}"
        )
        mapping_dubl.to_csv(duplicated_samples_name)
    else:
        print("WES<->microarray mapping is unique")


def verify_wes2microarray_mapping(
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
        "MicroARRAY mapping",
        results_file_prefix + ".MicroARRAY_mapping_non-unique.txt",
    )

    print("\n=== Mapping consistency checking ===")
    check_ids_consistency(ids_wes_data, ids_wes_mapping, "WES", results_file_prefix + ".WES")
    check_ids_consistency(
        ids_microarray_data, ids_microarray_mapping, "MicroARRAY", results_file_prefix + ".MicroARRAY"
    )

    mapping_marked_duplicated = mapping_mark_duplicates(mapping, wes_id_col, microarray_id_col)
    report_mapping_duplicated_stats(mapping_marked_duplicated, results_file_prefix + ".non-unique-pairs.csv")


def mapping_clean_missing(
    wes2microarray_mapping: pd.DataFrame,
    ids_microarray_data: IdList,
    microarray_id_col="microarray_id",
) -> pd.DataFrame:
    print("Removing from the mapping all microarray IDs not present in the real microarray data")
    ids_to_keep = wes2microarray_mapping[microarray_id_col].isin(ids_microarray_data)
    wes2array_clean_missing = wes2microarray_mapping.loc[ids_to_keep, :]
    ids_microarray_mapping_corrected_clean_missing = set(wes2array_clean_missing[microarray_id_col])
    print("Survived unique Array IDs: ", len(ids_microarray_mapping_corrected_clean_missing))
    return wes2array_clean_missing


# ======== The second part - validation of the gtcheck data consistency =======
def prepare_gtcheck(gtcheck: pd.DataFrame) -> pd.DataFrame:
    gtcheck_column_names = {
        i: s
        for i, s in enumerate(
            "DCv2 data_sample microarray_sample discordance average_logP n_sites N_matching_genotypes".split()
        )
    }
    gtcheck = gtcheck.rename(columns=gtcheck_column_names)
    gtcheck = gtcheck.drop("DCv2", axis=1)

    print("Checking mapping ID pairs uniqiness")

    if gtcheck.duplicated(subset=["data_sample", "microarray_sample"]).any():
        print("WARNING: Non-unique ID combinations found in the index. Duplicated items are dropped.")
        gtcheck = gtcheck.drop_duplicates(subset=["data_sample", "microarray_sample"], keep="first")
    else:
        print("ALL WES-microarray pairs are unique.")
    # Calculation score
    # Lowest score corresponds to the best match
    gtcheck["score"] = gtcheck.discordance / gtcheck.n_sites
    # Filling all non-calculated values with the "bad" score
    gtcheck["score"] = gtcheck["score"].replace({np.inf: 1000.0})
    return gtcheck


def prepare_mapping(mapping: pd.DataFrame, wes_id_col: str, microarray_id_col: str) -> pd.DataFrame:
    # Preparing mapping
    mapping = mapping.rename(columns={wes_id_col: "data_sample", microarray_id_col: "microarray_sample"})
    mapping = mapping[["data_sample", "microarray_sample"]]
    mapping = mapping_mark_duplicates(mapping, "data_sample", "microarray_sample")
    if mapping.duplicated_id_both.any():
        print(
            "WARNING - the mapping is double-way non unique!. Please double-check results of the mapping validation steps"
        )
    return mapping


def make_best_match(gtcheck: pd.DataFrame) -> pd.DataFrame:
    # Constructing table with the best score
    gtcheck_min_by_group_loc = gtcheck.groupby("data_sample").score.idxmin()
    gtcheck_best_match = gtcheck.loc[gtcheck_min_by_group_loc]
    gtcheck_best_match.set_index(keys=("data_sample"), inplace=True, drop=False)
    gtcheck_best_match.rename(
        columns={"microarray_sample": "best_match_microarray_sample", "score": "best_match_microarray_score"},
        inplace=True,
    )
    gtcheck_best_match.drop(["average_logP", "N_matching_genotypes"], axis=1, inplace=True)
    return gtcheck_best_match

    # The following set of functions separates the initial gtchech map to separate tables,
    # Each function returns the processed subtable with the marked validation result, and main table to process


def add_validation_result(df: pd.DataFrame, result: str, result_column_name: str = "validation_result"):
    """
    Adds a string result, to dataframe column result_column_name
    """
    df.loc[:, result_column_name] = df.loc[:, result_column_name] + result + ","


def set_validation_status(df: pd.DataFrame, status: bool, column_name="is_validation_passed") -> None:
    df.loc[:, column_name] = status


def filter_not_in_map(gtcheck_map: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Split gtcheck results to samples existing in map and not existing in map
    """
    mask_not_exist_in_map = gtcheck_map.microarray_sample.isnull()
    gtcheck_not_exist_in_map = gtcheck_map.loc[mask_not_exist_in_map]
    add_validation_result(gtcheck_not_exist_in_map, "best_match_not_exist_in_mapfile")
    gtcheck_exist_in_map = gtcheck_map.loc[~mask_not_exist_in_map]
    add_validation_result(gtcheck_exist_in_map, "best_match_exist_in_mapfile")
    return gtcheck_not_exist_in_map, gtcheck_exist_in_map


def filter_matched_map(gtcheck_map: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Splits gtcheck results to matching map and not matching mapping file
    """
    mask_no_pair_in_map = gtcheck_map.apply(
        lambda row: row.best_match_microarray_sample in row.microarray_sample, axis=1
    )

    gtcheck_matched_map = gtcheck_map.loc[mask_no_pair_in_map, :]
    add_validation_result(gtcheck_matched_map, "best_match_matched_mapfile")
    gtcheck_not_matched_map = gtcheck_map.loc[~mask_no_pair_in_map, :]
    add_validation_result(gtcheck_not_matched_map, "best_match_not_matched_mapfile")

    return gtcheck_matched_map, gtcheck_not_matched_map


def filter_by_score(
    df: pd.DataFrame, score_treshold: float, score_column_name: str = "score"
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """ """
    df_pass_score = df[df[score_column_name] <= score_treshold]
    add_validation_result(df_pass_score, "score_passed")
    df_fail_score = df[df[score_column_name] > score_treshold]
    add_validation_result(df_fail_score, "score_failed")
    return df_pass_score, df_fail_score


def mark_non_unique_matches(gtcheck_map: pd.DataFrame) -> pd.DataFrame:
    gtcheck_map_non_unique = gtcheck_map.loc[gtcheck_map.duplicated_id_any.eq(True), :]
    add_validation_result(gtcheck_map_non_unique, "mapfile_non_unique")
    gtcheck_map_unique = gtcheck_map.loc[gtcheck_map.duplicated_id_any.eq(False), :]
    add_validation_result(gtcheck_map_unique, "mapfile_unique")
    return pd.concat([gtcheck_map_unique, gtcheck_map_non_unique])


# TODO: make unit test
def validate_map_by_score(
    gtcheck_not_matched_map: pd.DataFrame, gtcheck: pd.DataFrame, mapping: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Analyses set of samples where best scores not found in mapping
    Finds the best score from all pairs existing in the mapping
    gtcheck - the original bcftools gtcheck results, containing all scores
    mapping - original mapping containing all pairs
    gtcheck_not_matched_map - samples that we're analysing
    """
    # Restoring the corresponding microarray sample by the mapping best match
    grouped_mapping = mapping.groupby(by="microarray_sample")["data_sample"].apply(list)
    mapping_dict = grouped_mapping.to_dict()
    wes_by_best_match = gtcheck_not_matched_map["best_match_microarray_sample"].map(mapping_dict)
    gtcheck_not_matched_map.loc[:, "wes_restored_by_best_match"] = wes_by_best_match.apply(
        lambda x: x if isinstance(x, list) else []
    )

    # Collecting for each WES possible mapping samples with their score
    gtcheck_not_matched_map_expanded = gtcheck_not_matched_map.explode("microarray_sample")

    gtcheck_score_dict = gtcheck.groupby(["data_sample", "microarray_sample"])["score"].first().to_dict()

    gtcheck_not_matched_map_expanded["score_from_mapped_microarray_sample"] = gtcheck_not_matched_map_expanded.apply(
        lambda row: gtcheck_score_dict.get((row.data_sample, row.microarray_sample), pd.NA), axis=1
    )

    have_no_gtscore = gtcheck_not_matched_map_expanded.score_from_mapped_microarray_sample.isnull()
    gtcheck_no_mapping_pairs = gtcheck_not_matched_map_expanded.loc[have_no_gtscore, :]
    add_validation_result(gtcheck_no_mapping_pairs, "no_mapfile_pairs_have_gtcheck")
    gtcheck_has_mapping_pairs = gtcheck_not_matched_map_expanded.loc[~have_no_gtscore, :]
    add_validation_result(gtcheck_has_mapping_pairs, "mapfile_pairs_have_gtcheck")
    return gtcheck_has_mapping_pairs, gtcheck_no_mapping_pairs


def gtcheck_validate(
    gtcheck: pd.DataFrame,
    mapping: pd.DataFrame,
    score_threshold: float = 0.2,
    wes_id_col: str = "data_sample",
    microarray_id_col: str = "microarray_sample",
) -> pd.DataFrame:
    gtcheck = prepare_gtcheck(gtcheck)
    gtcheck_best_match = make_best_match(gtcheck)

    mapping = prepare_mapping(mapping, wes_id_col, microarray_id_col)
    mapping_by_wes = mapping.groupby(by="data_sample").agg(
        {
            "microarray_sample": list,
            "duplicated_id_any": "any",
        }
    )
    print("Checkin for duplicated ids: ", mapping_by_wes.duplicated_id_any.value_counts())

    # Validation set of samples by mapping
    gtcheck_map = gtcheck_best_match.join(mapping_by_wes)
    with pd.option_context("future.no_silent_downcasting", True):
        gtcheck_map["duplicated_id_any"] = gtcheck_map.duplicated_id_any.fillna(False)

    # This column will contain the list of tags - text messages associated with each validation check
    gtcheck_map["validation_result"] = ""
    gtcheck_map["is_validation_passed"] = False

    # 1. WES IDs are not present in the map file - we can't do anything with it
    gtcheck_not_in_map, gtcheck_in_map = filter_not_in_map(gtcheck_map)
    set_validation_status(gtcheck_not_in_map, True)

    # 2. Wes-microarray pair exists in the map
    gtcheck_matched_map, gtcheck_not_matched_map = filter_matched_map(gtcheck_in_map)
    gtcheck_matched_map = mark_non_unique_matches(gtcheck_matched_map)

    # Filtering matched samples by score
    gtcheck_matched_map_pass_score, gtcheck_matched_map_fail_score = filter_by_score(
        gtcheck_matched_map, score_threshold, score_column_name="best_match_microarray_score"
    )
    set_validation_status(gtcheck_matched_map_pass_score, True)

    # 3. Recovering possible scores from map
    gtcheck_has_mapping_pairs, gtcheck_no_mapping_pairs = validate_map_by_score(
        gtcheck_not_matched_map, gtcheck, mapping
    )
    gtcheck_not_matched_pass_score, gtcheck_not_matched_fail_score = filter_by_score(
        gtcheck_has_mapping_pairs, score_threshold, score_column_name="score_from_mapped_microarray_sample"
    )
    set_validation_status(gtcheck_not_matched_pass_score, True)

    return pd.concat(
        [
            gtcheck_not_in_map,
            gtcheck_matched_map_pass_score,
            gtcheck_matched_map_fail_score,
            gtcheck_no_mapping_pairs,
            gtcheck_not_matched_pass_score,
            gtcheck_not_matched_fail_score,
        ]
    )


def print_validation_summary(df: pd.DataFrame) -> None:
    validation_result = df.validation_result.value_counts()
    print(validation_result)


def main() -> None:
    # = STEP SETUP = #
    config = parse_config()

    tmp_dir: str = config["general"]["tmp_dir"]

    # = STEP PARAMETERS = #
    dataset: str = config["general"]["dataset_name"]
    wes_id_col: str = config["step1"]["validate_gtcheck"]["wes_id_col"]
    microarray_id_col: str = config["step1"]["validate_gtcheck"]["microarray_id_col"]
    gtcheck_score_threshold: float = config["step1"]["validate_gtcheck"]["gtcheck_score_threshold"]

    # = STEP DEPENDENCIES = #
    mtpath: str = config["step1"]["validate_verifybamid"]["mt_metadata_annotated"]
    wes_microarray_mapping: Optional[str] = config["step1"]["validate_gtcheck"]["wes_microarray_mapping"]
    microarray_ids: str = config["step1"]["validate_gtcheck"]["microarray_ids"]
    gtcheck_report: str = config["step1"]["validate_gtcheck"]["gtcheck_report"]

    # = STEP OUTPUTS = #
    gtcheck_results_folder: str = config["step1"]["validate_gtcheck"]["gtcheck_results_folder"]
    gtcheck_results_prefix = os.path.join(gtcheck_results_folder, dataset)
    mt_gtcheck_validated: str = config["step1"]["validate_gtcheck"]["mt_gtcheck_validated"]

    # = STEP LOGIC = #
    _ = hail_utils.init_hl(tmp_dir)

    # Loading the list of WES samples from Hial matrixtable
    mt = hl.read_matrix_table(path_spark(mtpath))

    if wes_microarray_mapping is None:
        # No mapping file - nothing to check
        mt.write(path_spark(mt_gtcheck_validated), overwrite=True)
        return

    # === The first part - validation if the  WES-micro microarray mapping consistency ===

    Path(gtcheck_results_folder).mkdir(parents=True, exist_ok=True)

    ids_wes_data = mt.s.collect()

    # Loading the list of microarray samples from an external file
    with open(microarray_ids) as fmicroarray:
        ids_microarray_data = [s.strip() for s in fmicroarray.readlines()]

    report_id_stats(ids_wes_data, "WES data", f"{gtcheck_results_folder}/{dataset}.non-unique-ids.data.WES.txt")
    report_id_stats(
        ids_microarray_data, "MicroARRAY data", f"{gtcheck_results_folder}/{dataset}.non-unique-ids.data.microarray.txt"
    )

    # - Validating mapping table
    wes2microarray = pd.read_csv(wes_microarray_mapping, sep="\t")
    verify_wes2microarray_mapping(
        ids_wes_data,
        ids_microarray_data,
        wes2microarray,
        wes_id_col,
        microarray_id_col,
        results_file_prefix=gtcheck_results_prefix,
    )

    # - Correcting the mapping file - Removing microarray IDs not present in the real microarray data
    wes2microarray = mapping_clean_missing(wes2microarray, ids_microarray_data, microarray_id_col)

    # === The second part - validation of the gtcheck results ===
    gtcheck = pd.read_csv(gtcheck_report, sep="\t", header=None)

    validated = gtcheck_validate(
        gtcheck, wes2microarray, gtcheck_score_threshold, wes_id_col=wes_id_col, microarray_id_col=microarray_id_col
    )
    validated.to_csv(gtcheck_results_prefix + ".gtcheck_validation.all_samples.csv", index=False)

    # Printing summary for passed and failed validation file
    is_passed = validated.is_validation_passed
    validated_passed = validated.loc[is_passed, :]
    validated_failed = validated.loc[~is_passed, :]

    print(f"=== Passed validation: {len(validated_passed)} samples")
    print_validation_summary(validated_passed)
    print("")
    print(f"=== Failed validation: {len(validated_failed)} samples")
    if len(validated_failed) > 0:
        print_validation_summary(validated_failed)
        validated_failed.to_csv(gtcheck_results_prefix + ".gtcheck_validation.failed_samples.csv", index=False)

    # - Putting validation data in the main hail matrixtable
    # print(validated.dtypes) # TODO: the resulting pandas df containd many 'object' type fields instead of strings. Need to investigare and fix this.
    # Keeping only columns where Hail can correctly impute type
    validated = validated[
        [
            "data_sample",
            "best_match_microarray_sample",
            "best_match_microarray_score",
            "validation_result",
            "is_validation_passed",
        ]
    ]
    validated_table = hl.Table.from_pandas(validated, key="data_sample")
    mt = mt.annotate_cols(gtcheck_validation=validated_table[mt.s])

    mt.write(path_spark(mt_gtcheck_validated), overwrite=True)


if __name__ == "__main__":
    main()
