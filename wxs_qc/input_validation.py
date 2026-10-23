"""
Validation helpers for the pipeline input files
"""

import hail as hl

ALLOWED_SELF_REPORTED_SEX = ("male", "female", "undefined")


def validate_self_reported_sex(ht: hl.Table, field: str = "self_reported_sex") -> None:
    """
    Checks that the self-reported sex field contains only allowed values.

    Values are expected to be lowercased by the caller.
    Missing values are treated as invalid.

    Raises:
        ValueError: if any value is not in ALLOWED_SELF_REPORTED_SEX.
    """
    value_counts = ht.aggregate(hl.agg.counter(ht[field]))
    invalid_values = {value: count for value, count in value_counts.items() if value not in ALLOWED_SELF_REPORTED_SEX}
    if invalid_values:
        invalid_values_str = ", ".join(
            f"{'<missing>' if value is None else repr(value)}: {count} samples"
            for value, count in invalid_values.items()
        )
        raise ValueError(
            f"Invalid values in the '{field}' column: {invalid_values_str}. "
            f"Allowed values (case-insensitive): {', '.join(ALLOWED_SELF_REPORTED_SEX)}"
        )
