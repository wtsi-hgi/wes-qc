import os
import json
from unittest.mock import patch
import importlib

qc_step_4_1 = importlib.import_module("4-genotype_qc.1-compare_hard_filter_combinations")


def test_hard_filters_cache_filter_results(tmp_path):
    """Test the cache_filter_results decorator functionality"""

    # Create a mock function that would normally do expensive computations
    @qc_step_4_1.cache_filter_results
    def mock_process_filters(var_type: str, **kwargs) -> dict:
        # This function should only be called when cache miss
        return {"metric1": 100, "metric2": 200}

    # Setup test parameters
    filter_name = "test_filter"
    json_dump_folder = str(tmp_path)
    var_type = "snv"

    # First call - should compute and cache results
    with patch("builtins.print") as mock_print:
        result1 = mock_process_filters(filter_name=filter_name, json_dump_folder=json_dump_folder, var_type=var_type)

    # Verify first call behavior
    assert result1 == {"metric1": 100, "metric2": 200}
    cache_file = os.path.join(json_dump_folder, f"{var_type}_hardfilters_{filter_name}.json")
    assert os.path.exists(cache_file)

    # Verify the cached data
    with open(cache_file) as f:
        cached_data = json.load(f)
    assert cached_data == {filter_name: {"metric1": 100, "metric2": 200}}

    # Second call - should use cached results
    with patch("builtins.print") as mock_print:
        result2 = mock_process_filters(filter_name=filter_name, json_dump_folder=json_dump_folder, var_type=var_type)

    # Verify second call returns cached results
    assert result2 == {"metric1": 100, "metric2": 200}
    mock_print.assert_called_with(f"--- Checkpoint data loaded from file {cache_file}")


def test_hard_filters_cache_filter_results_different_params(tmp_path):
    """Test the cache_filter_results decorator with different parameters"""

    @qc_step_4_1.cache_filter_results
    def mock_process_filters(var_type: str, **kwargs) -> dict:
        return {"metric1": 100, "metric2": 200}

    # Call with different filter names
    filter_names = ["filter1", "filter2"]
    var_type = "snv"

    for filter_name in filter_names:
        _ = mock_process_filters(filter_name=filter_name, json_dump_folder=str(tmp_path), var_type=var_type)

        # Verify unique cache file created for each filter
        cache_file = os.path.join(tmp_path, f"{var_type}_hardfilters_{filter_name}.json")
        assert os.path.exists(cache_file)

        with open(cache_file) as f:
            cached_data = json.load(f)
        assert cached_data == {filter_name: {"metric1": 100, "metric2": 200}}
