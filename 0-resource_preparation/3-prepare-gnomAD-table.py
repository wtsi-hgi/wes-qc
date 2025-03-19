"""
This script prepares the gnomAD table,
keeping only data that are required for the VariantQC
"""

import hail as hl  # type: ignore

from utils.utils import parse_config
from utils.config import path_spark
from wes_qc import hail_utils


def prepare_gnomad_ht(ht: hl.Table) -> hl.Table:
    """
    Filter the Hail Table to keep only freq[0] and freq_meta[0].
    :param ht: Input Hail Table
    :return: Filtered Hail Table
    """
    ht_reduced = ht.select(freq=ht.freq[:2])
    ht_reduced = ht_reduced.select_globals(ht_reduced.freq_meta, ht_reduced.freq_index_dict)
    ht_reduced = ht_reduced.repartition(240)
    return ht_reduced


def main() -> None:
    # = STEP SETUP = #
    config = parse_config()

    tmp_dir = config["general"]["tmp_dir"]

    # = STEP PARAMETERS = #

    # = STEP DEPENDENCIES = #
    input_gnomad_htfile = config["step0"]["prepare_gnomad_ht"]["input_gnomad_htfile"]

    # = STEP OUTPUTS = #
    processed_gnomad_htfile = config["step0"]["prepare_gnomad_ht"]["processed_gnomad_htfile"]

    # = STEP LOGIC = #
    _ = hail_utils.init_hl(tmp_dir)
    input_gnomad_ht = hl.read_table(path_spark(input_gnomad_htfile))
    processed_gnomad_ht = prepare_gnomad_ht(input_gnomad_ht)
    processed_gnomad_ht.write(path_spark(processed_gnomad_htfile), overwrite=True)


if __name__ == "__main__":
    main()
