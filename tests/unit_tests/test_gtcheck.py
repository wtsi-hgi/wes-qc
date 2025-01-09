import importlib
import pandas as pd
from pytest import mark as m

qc_step_1_3 = importlib.import_module("1-import_data.3-validate-gtcheck")


@m.context("")
@m.it("")
def test_gtcheck_check_ids_consistency() -> None:
    pass


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
