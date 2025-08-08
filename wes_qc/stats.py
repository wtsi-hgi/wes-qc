"""
All functions used to calculate sample/variant statistics
"""

import hail as hl


def count_vars_per_cq(ht: hl.Table, cqs: list[str]) -> dict[str, int]:
    """
    Count variants per consequence type in one pass.
    Assumes ht.info.consequence is a string with consequences joined by '&'.
    """
    # Split consequences and explode into one row per consequence
    ht = ht.annotate(consequence_split=ht.info.consequence.split("&"))
    ht_exploded = ht.explode(ht.consequence_split)

    # Filter only consequences of interest
    ht_filtered = ht_exploded.filter(hl.literal(set(cqs)).contains(ht_exploded.consequence_split))

    # Group by consequence and count
    result_ht = ht_filtered.group_by(ht_filtered.consequence_split).aggregate(n=hl.agg.count())
    # Collect into Python dict
    result_dict = {row.consequence_split: row.n for row in result_ht.collect()}
    return result_dict
