# export to VCF after annotating with a range of hard filters - stringent filters only, remove other gts
from typing import Union

import hail as hl
import os.path
from utils.utils import parse_config, path_spark
from wes_qc import hail_utils


def export_vcfs(
    mtfile: str, filtered_vcf_dir: str, hard_filters: dict[str, dict[str, dict[str : Union[int, float]]]], model_id: str
):
    """
    Export VCFs annotated with a range of hard filters
    :param str mtfile: matrixtable file
    :param str filtered_vcf_dir: output directory for VCFs
    :param dict hard_filters: details of all sets of filters
    :param str model_id: random forest run hash used
    """
    mt = hl.read_matrix_table(path_spark(mtfile))

    # filter to remove rows where all variants fail the most stringent filters,
    mt = mt.filter_rows(mt.info.fraction_pass_stringent_filters > 0)
    mt = mt.filter_entries(mt.stringent_filters == "Pass")

    # drop unwanted fields
    mt = mt.drop(
        mt.a_index,
        mt.was_split,
        mt.variant_qc,
        mt.stringent_pass_count,
        mt.medium_pass_count,
        mt.relaxed_pass_count,
        mt.adj,
        mt.assigned_pop,
        mt.sum_AD,
    )
    mt = mt.annotate_rows(info=mt.info.drop("fraction_pass_medium_filters", "fraction_pass_relaxed_filters"))

    # info for header
    stringent_filters = (
        "SNPs: RF bin<="
        + str(hard_filters["snp"]["stringent"]["bin"])
        + " & DP>="
        + str(hard_filters["snp"]["stringent"]["dp"])
        + " & GQ>="
        + str(hard_filters["snp"]["stringent"]["gq"])
        + " & HetAB>="
        + str(hard_filters["snp"]["stringent"]["ab"])
        + " & Call_Rate>="
        + str(hard_filters["snp"]["stringent"]["call_rate"])
        + ", Indels: RF bin<="
        + str(hard_filters["indel"]["stringent"]["bin"])
        + " & DP>="
        + str(hard_filters["indel"]["stringent"]["dp"])
        + " & GQ>="
        + str(hard_filters["indel"]["stringent"]["gq"])
        + " & HetAB>="
        + str(hard_filters["indel"]["stringent"]["ab"])
        + " & Call_Rate>="
        + str(hard_filters["indel"]["stringent"]["call_rate"])
    )

    metadata = {
        "format": {
            "HetAB": {"Description": "Hetrozygous allele balance", "Number": "A", "Type": "Float"},
            "stringent_filters": {
                "Description": "Pass/fail stringent hard filters " + stringent_filters,
                "Number": "A",
                "Type": "String",
            },
        },
        "info": {
            "fraction_pass_stringent_filters": {
                "Description": "Fraction of genotypes which pass stringent hard filters " + stringent_filters,
                "Number": "A",
                "Type": "Float",
            },
            "rf_score": {
                "Description": "Variant QC random forest score, model id " + model_id,
                "Number": "A",
                "Type": "Float",
            },
            "rf_bin": {
                "Description": "Variant QC random forest bin, model id " + model_id,
                "Number": "A",
                "Type": "Integer",
            },
            "consequence": {"Description": "Most severe consequence from VEP104", "Number": "A", "Type": "String"},
            "gene": {
                "Description": "Gene affected by the most severe consequence from VEP104",
                "Number": "A",
                "Type": "String",
            },
            "hgnc_id": {
                "Description": "HGNC id of the gene affected by the most severe consequence from VEP104",
                "Number": "A",
                "Type": "String",
            },
        },
    }

    # export per chromosome
    chroms = [*range(1, 23), "X", "Y"]
    chromosomes = [f"chr{str(chr)}" for chr in chroms]
    for chromosome in chromosomes:
        print("Exporting " + chromosome)
        mt_chrom = mt.annotate_globals(chromosome=chromosome)
        mt_chrom = mt_chrom.filter_rows(mt_chrom.locus.contig == mt_chrom.chromosome)
        outfile = os.path.join(path_spark(filtered_vcf_dir), f"{chromosome}_stringent_filters.vcf.bgz")
        hl.export_vcf(mt_chrom, outfile, metadata=metadata)


def main():
    # = STEP SETUP = #
    config = parse_config()
    tmp_dir = config["general"]["tmp_dir"]

    # = STEP PARAMETERS = #
    model_id = config["general"]["rf_model_id"]
    hard_filters = config["step4"]["apply_hard_filters"]["hard_filters"]  # set during step 4.1a

    # = STEP DEPENDENCIES = #
    mtfile = config["step4"]["export_vcfs_b"]["mtfile"]

    # = STEP OUTPUTS = #
    filtered_vcf_dir = config["step4"]["export_vcfs_b"]["filtered_vcf_dir"]

    # = STEP LOGIC = #
    _ = hail_utils.init_hl(tmp_dir)
    export_vcfs(mtfile, filtered_vcf_dir, hard_filters, model_id)


if __name__ == "__main__":
    main()
