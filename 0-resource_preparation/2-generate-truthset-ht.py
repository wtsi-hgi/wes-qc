# generate truth set table for variant QC random forest
import hail as hl
from utils.utils import parse_config, path_spark
from wes_qc import hail_utils


def get_truth_ht(omni: str, mills: str, thousand_genomes: str, hapmap: str, **kwargs) -> hl.Table:
    """
    Generate truth set HT
    :param str omni: file containing kgp_omni 1000 Genomes intersection Onni 2.5M array
    :param str mills: file containing Mills & Devine indels
    :param str thousand_genomes: file containing high confidence sites in 1000 genonmes
    :param str hapmap: file containing hapmap
    :param str truth_ht_file: output file name
    """
    omni_ht = hl.read_table(path_spark(omni))
    mills_ht = hl.read_table(path_spark(mills))
    thousand_genomes_ht = hl.read_table(path_spark(thousand_genomes))
    hapmap_ht = hl.read_table(path_spark(hapmap))

    truth_ht = (
        hapmap_ht.select(hapmap=True)
        .join(omni_ht.select(omni=True), how="outer")
        .join(thousand_genomes_ht.select(kgp_phase1_hc=True), how="outer")
        .join(mills_ht.select(mills=True), how="outer")
        .repartition(200, shuffle=False)
    )  # TODO: I suppose this number should not be considered as a config parameter candidate, should it?
    return truth_ht


def main() -> None:
    # set up
    config = parse_config()
    tmp_dir = config["general"]["tmp_dir"]

    # = STEP PARAMETERS = #

    # = STEP DEPENDENCIES = #
    omni = config["step0"]["generate_truthset_ht"]["omni"]
    mills = config["step0"]["generate_truthset_ht"]["mills"]
    thousand_genomes = config["step0"]["generate_truthset_ht"]["thousand_genomes"]
    hapmap = config["step0"]["generate_truthset_ht"]["hapmap"]

    # = STEP OUTPUTS = #
    truth_ht_outfile = config["step0"]["generate_truthset_ht"]["truth_ht_outfile"]

    # = STEP LOGIC = #
    hail_utils.init_hl(tmp_dir)
    truth_ht = get_truth_ht(omni, mills, thousand_genomes, hapmap)
    truth_ht.write(path_spark(truth_ht_outfile), overwrite=True)


if __name__ == "__main__":
    main()
