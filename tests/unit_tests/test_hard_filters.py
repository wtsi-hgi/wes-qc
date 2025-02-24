import os
import json
from unittest.mock import patch
import importlib
import hail as hl

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


def test_hard_filters_calculate_normalized_mendel_errors_proband():
    """
    Test calculation of normalized Mendelian errors using synthetic data.
    Creates a dataset with controlled Mendelian inheritance patterns including errors.
    """
    # Create a synthetic pedigree with one trio
    trio = hl.Trio(s="child", pat_id="father", mat_id="mother", fam_id="fam1", is_female=True)
    ped = hl.Pedigree([trio])

    # Create synthetic genotype data
    # We'll create 10 variants:
    # - 8 variants with valid inheritance
    # - 2 variants with Mendelian errors
    n_variants = 10
    n_samples = 3  # father, mother, child

    # Create base MT structure
    mt = hl.balding_nichols_model(
        n_populations=1, n_samples=n_samples, n_variants=n_variants, reference_genome="GRCh38"
    )

    # Create sample mapping
    sample_ids = ["father", "mother", "child"]
    mt = mt.annotate_cols(s=hl.literal(sample_ids)[mt.sample_idx])
    mt = mt.key_cols_by(mt.s)

    # Add row index for variant access
    mt = mt.add_row_index()

    # Create specific genotypes for testing
    # Format: [father_gt, mother_gt, child_gt]
    # Valid inheritance patterns:
    valid_gts = [
        [hl.Call([0, 0]), hl.Call([0, 0]), hl.Call([0, 0])],  # All homozygous ref
        [hl.Call([0, 1]), hl.Call([0, 0]), hl.Call([0, 0])],  # Het father, hom ref mother and child
        [hl.Call([0, 0]), hl.Call([0, 1]), hl.Call([0, 1])],  # Hom ref father, het mother, het child
        [hl.Call([1, 1]), hl.Call([0, 1]), hl.Call([0, 1])],  # Hom alt father, het mother, het child
        [hl.Call([0, 1]), hl.Call([0, 1]), hl.Call([1, 1])],  # Het parents, hom alt child
        [hl.Call([0, 1]), hl.Call([0, 1]), hl.Call([0, 0])],  # Het parents, hom ref child
        [hl.Call([0, 1]), hl.Call([0, 1]), hl.Call([0, 1])],  # All het
        [hl.Call([1, 1]), hl.Call([1, 1]), hl.Call([1, 1])],  # All hom alt
    ]

    # Invalid inheritance patterns (Mendelian errors):
    invalid_gts = [
        [hl.Call([0, 0]), hl.Call([0, 0]), hl.Call([0, 1])],  # Hom ref parents, het child
        [hl.Call([0, 0]), hl.Call([0, 0]), hl.Call([1, 1])],  # Hom ref parents, hom alt child
    ]

    # Combine all genotypes
    all_gts = valid_gts + invalid_gts

    # Create a list of genotypes per sample
    father_gts = [gt[0] for gt in all_gts]
    mother_gts = [gt[1] for gt in all_gts]
    child_gts = [gt[2] for gt in all_gts]

    # Create an annotation expression for the genotypes
    mt = mt.annotate_entries(
        GT=hl.case()
        .when(mt.s == "father", hl.literal(father_gts)[hl.int32(mt.row_idx)])
        .when(mt.s == "mother", hl.literal(mother_gts)[hl.int32(mt.row_idx)])
        .when(mt.s == "child", hl.literal(child_gts)[hl.int32(mt.row_idx)])
        .or_missing()
    )

    # Calculate normalized Mendelian errors
    mean_error_rate, std_error_rate = qc_step_4_1.calculate_normalized_mendel_errors_proband(mt, ped)

    # Expected error rate: 2 errors out of 10 variants = 0.2
    assert abs(mean_error_rate - 0.2) < 1e-6, f"Expected error rate 0.2, got {mean_error_rate}"
    # Since we only have one trio, std dev should be 0
    assert abs(std_error_rate) < 1e-6, f"Expected std dev 0, got {std_error_rate}"
