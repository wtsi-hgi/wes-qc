import importlib
import pandas as pd
from pytest import mark as m

qc_step_1_3 = importlib.import_module("1-import_data.3-validate-gtcheck")


@m.context("Given lists of IDs from data and mapping")
@m.it("Should correctly identify missing IDs and handle duplicates")
def test_gtcheck_check_ids_consistency(tmp_path) -> None:
    # Test data
    data_ids = ["sample1", "sample2", "sample3", "sample4", "sample4"]  # Note: sample4 is duplicated
    map_ids = ["sample2", "sample3", "sample5", "sample5"]  # Note: sample5 is duplicated

    # Expected results
    expected_data_no_map = {"sample1", "sample4"}  # IDs in data but not in map
    expected_map_no_data = {"sample5"}  # IDs in map but not in data

    # Call the function
    data_no_map, map_no_data = qc_step_1_3.check_ids_consistency(
        data_ids=data_ids, map_ids=map_ids, caption="test", dump_prefix=str(tmp_path / "test")
    )

    # Check results
    assert data_no_map == expected_data_no_map
    assert map_no_data == expected_map_no_data

    # Check if dump files were created when mismatches found
    assert (tmp_path / "test.ids_data_not_in_map.txt").exists()
    assert (tmp_path / "test.ids_map_not_in_data.txt").exists()

    # Test case with perfect match
    perfect_data = ["sample1", "sample2"]
    perfect_map = ["sample1", "sample2"]

    data_no_map, map_no_data = qc_step_1_3.check_ids_consistency(
        data_ids=perfect_data, map_ids=perfect_map, caption="perfect", dump_prefix=str(tmp_path / "perfect")
    )

    # Check results for perfect match
    assert len(data_no_map) == 0
    assert len(map_no_data) == 0

    # Check that no dump files were created for perfect match
    assert not (tmp_path / "perfect.ids_data_not_in_map.txt").exists()
    assert not (tmp_path / "perfect.ids_map_not_in_data.txt").exists()


def test_gtcheck_mapping_mark_duplicates() -> None:
    # Construct a DataFrame with some duplicates
    mapping = pd.DataFrame({"wes_id": ["W1", "W2", "W3", "W1", "W4"], "microarray_id": ["A1", "A2", "A3", "A4", "A1"]})

    # Call the function
    result = qc_step_1_3.mapping_mark_duplicates(mapping)

    print("\n", result["duplicated_wes"])

    # Check correctness of `duplicated_*` columns
    assert (result["duplicated_wes"] == [True, False, False, True, False]).all()
    assert (result["duplicated_microarray"] == [True, False, False, False, True]).all()
    assert (result["duplicated_id_any"] == [True, False, False, True, True]).all()
    assert (result["duplicated_id_both"] == [True, False, False, False, False]).all()

    # Check if result still has the original columns
    assert "wes_id" in result
    assert "microarray_id" in result

    # Check if function does not modify the input DataFrame
    assert "duplicated_wes" not in mapping
    assert "duplicated_microarray" not in mapping
    assert "duplicated_id_any" not in mapping
    assert "duplicated_id_both" not in mapping

    # Check if dataframe 'mapping' still has the same 'wes_id' and 'array_id' data
    pd.testing.assert_series_equal(mapping["wes_id"], pd.Series(["W1", "W2", "W3", "W1", "W4"], name="wes_id"))
    pd.testing.assert_series_equal(
        mapping["microarray_id"], pd.Series(["A1", "A2", "A3", "A4", "A1"], name="microarray_id")
    )


@m.context("Given a list of sample IDs")
@m.it("Should correctly report stats for unique and non-unique IDs")
def test_gtcheck_report_id_stats(tmp_path, capsys) -> None:
    # Test case 1: List with unique IDs
    unique_ids = ["sample1", "sample2", "sample3"]
    stats_file = tmp_path / "unique_stats.txt"

    qc_step_1_3.report_id_stats(unique_ids, "unique test", stats_file)

    # Capture printed output
    captured = capsys.readouterr()
    expected_output = (
        "\n=== unique test samples validation ===\n"
        "Total 3 IDs in unique test data\n"
        "3 unique IDs in unique test\n"
        "All IDs in unique test are unique\n"
    )
    assert captured.out == expected_output
    # Check that stats file wasn't created (no non-unique IDs)
    assert not stats_file.exists()

    # Test case 2: List with non-unique IDs
    non_unique_ids = ["sample1", "sample2", "sample2", "sample3", "sample3"]
    stats_file = tmp_path / "non_unique_stats.txt"

    qc_step_1_3.report_id_stats(non_unique_ids, "non-unique test", stats_file)

    # Capture printed output
    captured = capsys.readouterr()
    expected_output = (
        "\n=== non-unique test samples validation ===\n"
        "Total 5 IDs in non-unique test data\n"
        "1 unique IDs in non-unique test\n"
        "2 non-unique IDs in non-unique test\n"
        f"Non-unique samples dumped to {stats_file}\n"
    )
    assert captured.out == expected_output

    # Check that stats file was created and contains correct non-unique IDs
    assert stats_file.exists()
    with open(stats_file) as f:
        non_unique_content = f.read().strip().split("\n")
    assert set(non_unique_content) == {"sample2", "sample3"}


@m.context("Given a mapping DataFrame with microarray IDs")
@m.it("Should remove entries where microarray ID is not in the provided data list")
def test_gtcheck_mapping_clean_missing() -> None:
    # Test data
    mapping = pd.DataFrame({"wes_id": ["W1", "W2", "W3", "W4", "W5"], "microarray_id": ["A1", "A2", "A3", "A4", "A5"]})

    # Only A1, A3, and A5 exist in actual data
    microarray_data_ids = ["A1", "A3", "A5"]

    # Call the function
    result = qc_step_1_3.mapping_clean_missing(
        wes2microarray_mapping=mapping, ids_microarray_data=microarray_data_ids, microarray_id_col="microarray_id"
    )

    # Check that only rows with existing microarray IDs remain
    expected = pd.DataFrame({"wes_id": ["W1", "W3", "W5"], "microarray_id": ["A1", "A3", "A5"]}, index=[0, 2, 4])

    pd.testing.assert_frame_equal(result, expected)

    # Test with empty microarray data list
    empty_result = qc_step_1_3.mapping_clean_missing(
        wes2microarray_mapping=mapping, ids_microarray_data=[], microarray_id_col="microarray_id"
    )
    assert len(empty_result) == 0

    # Test with custom column name
    mapping_custom = pd.DataFrame({"wes_id": ["W1", "W2"], "custom_array_col": ["A1", "A2"]})
    custom_result = qc_step_1_3.mapping_clean_missing(
        wes2microarray_mapping=mapping_custom, ids_microarray_data=["A1"], microarray_id_col="custom_array_col"
    )
    assert len(custom_result) == 1
    assert custom_result.iloc[0]["custom_array_col"] == "A1"


@m.context("Given a gtcheck DataFrame from bcftools")
@m.it("Should prepare the DataFrame by renaming columns, checking uniqueness, and calculating scores")
def test_gtcheck_prepare_gtcheck(capsys) -> None:
    # Test data
    gtcheck = pd.DataFrame(
        {
            0: ["DCv2"] * 4,
            1: ["sample1", "sample2", "sample2", "sample3"],  # data_sample
            2: ["array1", "array2", "array2", "array3"],  # microarray_sample
            3: [0.01, 0.02, 0.02, 0.009],  # discordance
            4: [-10.5, -9.8, -9.8, -8.5],  # average_logP
            5: [100, 95, 95, 0],  # n_sites
            6: [98, 93, 93, 88],  # N_matching_genotypes
        }
    )

    # Call the function
    result = qc_step_1_3.prepare_gtcheck(gtcheck)

    # Capture printed output
    captured = capsys.readouterr()
    assert "Checking mapping ID pairs uniqiness" in captured.out
    assert "WARNING: Non-unique ID combinations found in the index" in captured.out

    # Check column renaming
    expected_columns = [
        "data_sample",
        "microarray_sample",
        "discordance",
        "average_logP",
        "n_sites",
        "N_matching_genotypes",
        "score",
    ]
    assert list(result.columns) == expected_columns

    # Check duplicates were removed
    assert len(result) == 3  # One duplicate pair should be removed
    assert not result.duplicated(subset=["data_sample", "microarray_sample"]).any()

    # Check score calculation
    expected_scores = [0.0001, 0.00021052631578947367, 1000]  # discordance/n_sites, 1 for NA
    pd.testing.assert_series_equal(
        result["score"].round(10),
        pd.Series(expected_scores, name="score", index=[0, 1, 3]).round(10),
        check_names=False,
    )

    # Test case with no duplicates
    unique_gtcheck = pd.DataFrame(
        {
            0: ["DCv2"] * 2,
            1: ["sample1", "sample2"],
            2: ["array1", "array2"],
            3: [0.01, 0.02],
            4: [-10.5, -9.8],
            5: [100, 95],
            6: [98, 93],
        }
    )

    _ = qc_step_1_3.prepare_gtcheck(unique_gtcheck)
    captured = capsys.readouterr()
    assert "ALL WES-microarray pairs are unique" in captured.out


@m.context("Given a mapping DataFrame")
@m.it("Should prepare the mapping by renaming columns and marking duplicates")
def test_gtcheck_prepare_mapping(capsys) -> None:
    # Test data with some duplicates
    mapping = pd.DataFrame(
        {
            "wes_col": ["W1", "W2", "W3", "W1", "W4"],
            "microarray_col": ["A1", "A2", "A3", "A4", "A1"],
            "extra_col": ["X1", "X2", "X3", "X4", "X5"],  # Should be dropped
        }
    )

    # Call the function
    result = qc_step_1_3.prepare_mapping(mapping=mapping, wes_id_col="wes_col", microarray_id_col="microarray_col")

    # Capture printed output
    captured = capsys.readouterr()
    assert "WARNING - the mapping is double-way non unique!" in captured.out

    # Check column renaming and selection
    expected_columns = [
        "data_sample",
        "microarray_sample",
        "duplicated_wes",
        "duplicated_microarray",
        "duplicated_id_any",
        "duplicated_id_both",
    ]
    assert list(result.columns) == expected_columns

    # Check if data was correctly renamed
    pd.testing.assert_series_equal(result["data_sample"], pd.Series(["W1", "W2", "W3", "W1", "W4"], name="data_sample"))
    pd.testing.assert_series_equal(
        result["microarray_sample"], pd.Series(["A1", "A2", "A3", "A4", "A1"], name="microarray_sample")
    )

    # Check duplicate marking
    assert (result["duplicated_wes"] == [True, False, False, True, False]).all()
    assert (result["duplicated_microarray"] == [True, False, False, False, True]).all()
    assert (result["duplicated_id_any"] == [True, False, False, True, True]).all()
    assert (result["duplicated_id_both"] == [True, False, False, False, False]).all()

    # Test case with no duplicates
    unique_mapping = pd.DataFrame({"wes_col": ["W1", "W2"], "microarray_col": ["A1", "A2"], "extra_col": ["X1", "X2"]})

    unique_result = qc_step_1_3.prepare_mapping(
        mapping=unique_mapping,
        wes_id_col="wes_col",
        microarray_id_col="microarray_col",  # Fixed to use consistent column name
    )

    # Capture printed output for unique case
    captured = capsys.readouterr()
    assert "WARNING" not in captured.out  # No warning should be printed for unique mappings

    # Check that no duplicates were marked
    assert not unique_result["duplicated_wes"].any()
    assert not unique_result["duplicated_microarray"].any()
    assert not unique_result["duplicated_id_any"].any()
    assert not unique_result["duplicated_id_both"].any()


@m.context("Given a gtcheck DataFrame with multiple matches per sample")
@m.it("Should select the best match based on lowest score for each data sample")
def test_gtcheck_make_best_match() -> None:
    # Test data with multiple matches per sample
    gtcheck = pd.DataFrame(
        {
            "data_sample": ["sample1", "sample1", "sample2", "sample2", "sample3"],
            "microarray_sample": ["array1", "array2", "array3", "array4", "array5"],
            "score": [0.1, 0.2, 0.05, 0.15, 0.3],
            "average_logP": [-10, -9, -8, -7, -6],
            "N_matching_genotypes": [95, 90, 98, 92, 85],
            "n_sites": [100, 100, 100, 100, 100],
            "discordance": [0.1, 0.2, 0.05, 0.15, 0.3],
        }
    )

    # Call the function
    result = qc_step_1_3.make_best_match(gtcheck)

    # Check that we got the best match for each sample
    assert len(result) == 3  # Should have one row per unique data_sample
    assert list(result.index) == ["sample1", "sample2", "sample3"]

    # Check that we selected the lowest scores
    assert result.loc["sample1", "best_match_microarray_score"] == 0.1  # array1 was better match
    assert result.loc["sample2", "best_match_microarray_score"] == 0.05  # array3 was better match
    assert result.loc["sample3", "best_match_microarray_score"] == 0.3  # only one match

    # Check that we kept the right columns
    expected_columns = {
        "data_sample",
        "best_match_microarray_sample",
        "best_match_microarray_score",
        "n_sites",
        "discordance",
    }
    assert set(result.columns) == expected_columns

    # Check that we dropped the right columns
    assert "average_logP" not in result.columns
    assert "N_matching_genotypes" not in result.columns


@m.context("Given a gtcheck DataFrame with mapping results")
@m.it("Should separate samples into those existing in map and those not existing")
def test_gtcheck_filter_not_in_map() -> None:
    # Test data with some samples having no mapping (null microarray_sample)
    gtcheck_map = pd.DataFrame(
        {
            "data_sample": ["sample1", "sample2", "sample3", "sample4"],
            "best_match_microarray_sample": ["array1", None, "array3", "array4"],
            "best_match_microarray_score": [0.1, 0.2, 0.15, 0.3],
            "microarray_sample": [["array1", "array5"], None, ["array3"], ["array6"]],
            "validation_result": ["", "", "", ""],
        }
    )

    # Call the function
    not_in_map, exist_in_map = qc_step_1_3.filter_not_in_map(gtcheck_map)

    # Check separation of samples
    assert len(not_in_map) == 1  # Only sample2 has no mapping
    assert len(exist_in_map) == 3  # The rest have mappings

    # Check that samples were correctly separated
    assert not_in_map.iloc[0]["data_sample"] == "sample2"
    assert set(exist_in_map["data_sample"]) == {"sample1", "sample3", "sample4"}

    # Check validation results were properly added
    assert "best_match_not_exist_in_mapfile" in not_in_map.iloc[0]["validation_result"]
    assert all("best_match_exist_in_mapfile" in result for result in exist_in_map["validation_result"])

    # Test case with all samples having mapping
    all_mapped = gtcheck_map.copy()
    all_mapped["microarray_sample"] = [["array1"], ["array2"], ["array3"], ["array4"]]
    all_mapped["validation_result"] = ["", "", "", ""]

    not_in_map, exist_in_map = qc_step_1_3.filter_not_in_map(all_mapped)
    assert len(not_in_map) == 0
    assert len(exist_in_map) == 4

    # Test case with no samples having mapping
    none_mapped = gtcheck_map.copy()
    none_mapped["microarray_sample"] = [None, None, None, None]
    none_mapped["validation_result"] = ["", "", "", ""]

    not_in_map, exist_in_map = qc_step_1_3.filter_not_in_map(none_mapped)
    assert len(not_in_map) == 4
    assert len(exist_in_map) == 0


@m.context("Given a gtcheck DataFrame with mapping results")
@m.it("Should separate samples into those matching and not matching the mapping file")
def test_gtcheck_filter_matched_map() -> None:
    # Test data with some samples matching their mapping and others not
    gtcheck_map = pd.DataFrame(
        {
            "data_sample": ["sample1", "sample2", "sample3", "sample4"],
            "best_match_microarray_sample": ["array1", "array2", "array3", "array5"],
            "best_match_microarray_score": [0.1, 0.2, 0.15, 0.3],
            "microarray_sample": [
                ["array1", "array5"],  # Matches best_match
                ["array3", "array4"],  # Doesn't match best_match
                ["array3"],  # Matches best_match
                ["array4"],  # Doesn't match best_match
            ],
            "validation_result": ["", "", "", ""],
        }
    )

    # Call the function
    matched_map, not_matched_map = qc_step_1_3.filter_matched_map(gtcheck_map)

    # Check separation of samples
    assert len(matched_map) == 2  # sample1 and sample3 match their best_match
    assert len(not_matched_map) == 2  # sample2 and sample4 don't match

    # Check that samples were correctly separated
    assert set(matched_map["data_sample"]) == {"sample1", "sample3"}
    assert set(not_matched_map["data_sample"]) == {"sample2", "sample4"}

    # Check validation results were properly added
    assert all("best_match_matched_mapfile" in result for result in matched_map["validation_result"])
    assert all("best_match_not_matched_mapfile" in result for result in not_matched_map["validation_result"])

    # Test case with all samples matching their mapping
    all_matched = gtcheck_map.copy()
    all_matched["microarray_sample"] = [["array1"], ["array2"], ["array3"], ["array5"]]
    all_matched["validation_result"] = ["", "", "", ""]

    matched_map, not_matched_map = qc_step_1_3.filter_matched_map(all_matched)
    assert len(matched_map) == 4
    assert len(not_matched_map) == 0

    # Test case with no samples matching their mapping
    none_matched = gtcheck_map.copy()
    none_matched["microarray_sample"] = [["arrayX"], ["arrayY"], ["arrayZ"], ["arrayW"]]
    none_matched["validation_result"] = ["", "", "", ""]

    matched_map, not_matched_map = qc_step_1_3.filter_matched_map(none_matched)
    assert len(matched_map) == 0
    assert len(not_matched_map) == 4


@m.context("Given a DataFrame with scores")
@m.it("Should separate records based on score threshold and add validation results")
def test_gtcheck_filter_by_score() -> None:
    # Test data with various scores
    df = pd.DataFrame(
        {
            "data_sample": ["sample1", "sample2", "sample3", "sample4", "sample5"],
            "score": [0.1, 0.2, 0.3, 0.15, 0.25],
            "validation_result": ["", "", "", "", ""],
        }
    )

    # Call the function with threshold 0.2
    passed_score, failed_score = qc_step_1_3.filter_by_score(df=df, score_treshold=0.2, score_column_name="score")

    # Check separation of samples
    assert len(passed_score) == 3  # sample1, sample2, sample4 should pass
    assert len(failed_score) == 2  # sample3, sample5 should fail

    # Check that samples were correctly separated based on score
    assert set(passed_score["data_sample"]) == {"sample1", "sample2", "sample4"}
    assert set(failed_score["data_sample"]) == {"sample3", "sample5"}

    # Check validation results were properly added
    assert all("score_passed" in result for result in passed_score["validation_result"])
    assert all("score_failed" in result for result in failed_score["validation_result"])

    # Test case with all samples passing threshold
    all_pass = df.copy()
    all_pass["score"] = [0.1, 0.1, 0.1, 0.1, 0.1]
    all_pass["validation_result"] = ["", "", "", "", ""]

    passed_score, failed_score = qc_step_1_3.filter_by_score(df=all_pass, score_treshold=0.2, score_column_name="score")
    assert len(passed_score) == 5
    assert len(failed_score) == 0

    # Test case with all samples failing threshold
    all_fail = df.copy()
    all_fail["score"] = [0.3, 0.3, 0.3, 0.3, 0.3]
    all_fail["validation_result"] = ["", "", "", "", ""]

    passed_score, failed_score = qc_step_1_3.filter_by_score(df=all_fail, score_treshold=0.2, score_column_name="score")
    assert len(passed_score) == 0
    assert len(failed_score) == 5

    # Test with custom score column name
    df_custom = pd.DataFrame(
        {"data_sample": ["sample1", "sample2"], "custom_score": [0.1, 0.3], "validation_result": ["", ""]}
    )

    passed_score, failed_score = qc_step_1_3.filter_by_score(
        df=df_custom, score_treshold=0.2, score_column_name="custom_score"
    )
    assert len(passed_score) == 1
    assert len(failed_score) == 1
    assert passed_score.iloc[0]["data_sample"] == "sample1"
    assert failed_score.iloc[0]["data_sample"] == "sample2"


@m.context("Given a gtcheck mapping DataFrame")
@m.it("Should correctly mark unique and non-unique matches and add validation results")
def test_gtcheck_mark_non_unique_matches() -> None:
    # Test data with mix of unique and non-unique mappings
    gtcheck_map = pd.DataFrame(
        {
            "data_sample": ["sample1", "sample2", "sample3", "sample4", "sample5"],
            "microarray_sample": ["array1", "array2", "array3", "array1", "array5"],
            "duplicated_id_any": [True, False, False, True, False],
            "validation_result": ["", "", "", "", ""],
        }
    )

    # Call the function
    result = qc_step_1_3.mark_non_unique_matches(gtcheck_map)

    # Check that all original rows are preserved
    assert len(result) == len(gtcheck_map)

    # Check that samples were correctly marked
    non_unique_samples = result[result["duplicated_id_any"]]
    unique_samples = result[~result["duplicated_id_any"]]

    # Explicit checks for sample names in non-unique group
    assert set(non_unique_samples["data_sample"]) == {"sample1", "sample4"}
    assert set(non_unique_samples["microarray_sample"]) == {"array1"}

    # Explicit checks for sample names in unique group
    assert set(unique_samples["data_sample"]) == {"sample2", "sample3", "sample5"}
    assert set(unique_samples["microarray_sample"]) == {"array2", "array3", "array5"}

    # Check counts
    assert len(non_unique_samples) == 2  # sample1 and sample4
    assert len(unique_samples) == 3  # sample2, sample3, and sample5

    # Check validation results were properly added
    assert all("mapfile_non_unique" in result for result in non_unique_samples["validation_result"])
    assert all("mapfile_unique" in result for result in unique_samples["validation_result"])

    # Test case with all unique mappings
    all_unique = gtcheck_map.copy()
    all_unique["duplicated_id_any"] = False
    all_unique["validation_result"] = ["", "", "", "", ""]

    result_all_unique = qc_step_1_3.mark_non_unique_matches(all_unique)
    # Explicit check for all sample names in unique group
    assert set(result_all_unique["data_sample"]) == {"sample1", "sample2", "sample3", "sample4", "sample5"}
    assert all("mapfile_unique" in result for result in result_all_unique["validation_result"])
    assert not any("mapfile_non_unique" in result for result in result_all_unique["validation_result"])

    # Test case with all non-unique mappings
    all_non_unique = gtcheck_map.copy()
    all_non_unique["duplicated_id_any"] = True
    all_non_unique["validation_result"] = ["", "", "", "", ""]

    result_all_non_unique = qc_step_1_3.mark_non_unique_matches(all_non_unique)
    # Explicit check for all sample names in non-unique group
    assert set(result_all_non_unique["data_sample"]) == {"sample1", "sample2", "sample3", "sample4", "sample5"}
    assert all("mapfile_non_unique" in result for result in result_all_non_unique["validation_result"])
    assert not any("mapfile_unique" in result for result in result_all_non_unique["validation_result"])


@m.context("Given gtcheck results with unmatched mapping")
@m.it("Should validate samples by comparing scores from mapping pairs")
def test_gtcheck_validate_map_by_score() -> None:
    # Test data setup
    gtcheck_not_matched_map = pd.DataFrame(
        {
            "data_sample": ["sample1", "sample2"],
            "best_match_microarray_sample": ["array2", "array4"],
            "microarray_sample": [["array1", "array3"], ["array3", "array5"]],
            "validation_result": ["", ""],
        }
    )

    # Original gtcheck results containing all scores
    gtcheck = pd.DataFrame(
        {
            "data_sample": ["sample1", "sample1", "sample1", "sample2", "sample2", "sample2"],
            "microarray_sample": ["array1", "array2", "array3", "array3", "array4", "array5"],
            "score": [0.15, 0.1, 0.2, 0.3, 0.25, 0.4],
        }
    )

    # Original mapping
    mapping = pd.DataFrame(
        {
            "data_sample": ["sample1", "sample1", "sample2", "sample2", "sample4"],
            "microarray_sample": ["array1", "array3", "array3", "array5", "array4"],
        }
    )

    # Call the function
    has_mapping_pairs, no_mapping_pairs = qc_step_1_3.validate_map_by_score(
        gtcheck_not_matched_map=gtcheck_not_matched_map, gtcheck=gtcheck, mapping=mapping
    )

    # Check results for samples with mapping pairs
    assert len(has_mapping_pairs) == 4  # sample1-array1, sample1-array3, sample2-array3, sample2-array5
    # Check results for samples without mapping pairs
    assert len(no_mapping_pairs) == 0  # All pairs had scores in this test case

    assert "mapfile_pairs_have_gtcheck" in has_mapping_pairs["validation_result"].iloc[0]

    # Verify scores were correctly matched
    scores = has_mapping_pairs["score_from_mapped_microarray_sample"].dropna().tolist()
    expected_scores = [0.15, 0.2, 0.3, 0.4]  # Scores for the mapped pairs
    assert scores == expected_scores

    # Test case with missing gtcheck scores
    gtcheck_missing = gtcheck.copy()
    gtcheck_missing = gtcheck_missing[gtcheck_missing.microarray_sample != "array1"]

    has_mapping_pairs, no_mapping_pairs = qc_step_1_3.validate_map_by_score(
        gtcheck_not_matched_map=gtcheck_not_matched_map, gtcheck=gtcheck_missing, mapping=mapping
    )

    # Verify handling of missing scores
    assert len(has_mapping_pairs) == 3  # Only pairs with scores
    assert len(no_mapping_pairs) == 1  # One pair missing score
    assert "no_mapfile_pairs_have_gtcheck" in no_mapping_pairs["validation_result"].iloc[0]
