# export to VCF after annotating with a range of hard filters
import datetime
import os
from typing import Any

import hail as hl
from utils.utils import parse_config
from wes_qc import hail_utils

fail_string = "FAIL"
outside_bait_string = "OUTSIDE_BAIT"
pass_medium_string = "PASS_MEDIUM"
pass_stringent_string = "PASS_STRINGENT"


def export_vcfs(
    mtfile: str, filtered_vcf_dir: str, hard_filters: dict[str, Any], run_hash: str, export_trace: bool = True
) -> None:
    """
    Export VCFs annotated with a range of hard filters
    :param str mtfile: matrixtable file
    :param str filtered_vcf_dir: output directory for VCFs
    :param dict hard_filters: details of all sets of filters
    :param str run_hash: random forest run hash used
    """
    mt = hl.read_matrix_table(mtfile)
    print(f"Exporting VCfs from matrixtable {mtfile}")

    # filter to remove rows where all variants fail the most relaxed filters
    mt = mt.filter_rows(mt.info.fraction_pass_relaxed_filters > 0)

    # drop unwanted fields
    mt = mt.annotate_rows(info=mt.info.drop("AC", "AN", "AF"))
    for filter_level in ("relaxed", "medium", "stringent"):
        for metric in ("AN", "AC", "AC_Hom", "AC_Het"):
            filter_metric = f"{filter_level}_{metric}"
            mt = mt.annotate_rows(info=mt.info.annotate(**{filter_metric: mt[filter_metric]}))
            mt = mt.drop(mt[filter_metric])

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
        + ", Indels: RF bin<="
        + str(hard_filters["indel"]["stringent"]["bin"])
        + " & DP>="
        + str(hard_filters["indel"]["stringent"]["dp"])
        + " & GQ>="
        + str(hard_filters["indel"]["stringent"]["gq"])
        + " & HetAB>="
        + str(hard_filters["indel"]["stringent"]["ab"])
        + " & Call Rate>="
        + str(hard_filters["missingness"])
    )

    medium_filters = (
        "SNPs: RF bin<="
        + str(hard_filters["snp"]["medium"]["bin"])
        + " & DP>="
        + str(hard_filters["snp"]["medium"]["dp"])
        + " & GQ>="
        + str(hard_filters["snp"]["medium"]["gq"])
        + " & HetAB>="
        + str(hard_filters["snp"]["medium"]["ab"])
        + ", Indels: RF bin<="
        + str(hard_filters["indel"]["medium"]["bin"])
        + " & DP>="
        + str(hard_filters["indel"]["medium"]["dp"])
        + " & GQ>="
        + str(hard_filters["indel"]["medium"]["gq"])
        + " & HetAB>="
        + str(hard_filters["indel"]["medium"]["ab"])
        + " & Call Rate>="
        + str(hard_filters["missingness"])
    )

    relaxed_filters = (
        "SNPs: RF bin<="
        + str(hard_filters["snp"]["relaxed"]["bin"])
        + " & DP>="
        + str(hard_filters["snp"]["relaxed"]["dp"])
        + " & GQ>="
        + str(hard_filters["snp"]["relaxed"]["gq"])
        + " & HetAB>="
        + str(hard_filters["snp"]["relaxed"]["ab"])
        + ", Indels: RF bin<="
        + str(hard_filters["indel"]["relaxed"]["bin"])
        + " & DP>="
        + str(hard_filters["indel"]["relaxed"]["dp"])
        + " & GQ>="
        + str(hard_filters["indel"]["relaxed"]["gq"])
        + " & HetAB>="
        + str(hard_filters["indel"]["relaxed"]["ab"])
        + " & Call Rate>="
        + str(hard_filters["missingness"])
    )

    metadata = {
        "info": {
            "stringent_AN": {
                "Description": "Total number of alleles in called genotypes",
                "Number": "1",
                "Type": "Integer",
            },
            "stringent_AC": {"Description": "Allele count in genotypes", "Number": "A", "Type": "Integer"},
            "stringent_AC_Hom": {
                "Description": "Allele counts in homozygous genotypes",
                "Number": "A",
                "Type": "Integer",
            },
            "stringent_AC_Het": {
                "Description": "Allele counts in heterozygous genotypes",
                "Number": "A",
                "Type": "Integer",
            },
            "medium_AN": {
                "Description": "Total number of alleles in called genotypes",
                "Number": "1",
                "Type": "Integer",
            },
            "medium_AC": {"Description": "Allele count in genotypes", "Number": "A", "Type": "Integer"},
            "medium_AC_Hom": {"Description": "Allele counts in homozygous genotypes", "Number": "A", "Type": "Integer"},
            "medium_AC_Het": {
                "Description": "Allele counts in heterozygous genotypes",
                "Number": "A",
                "Type": "Integer",
            },
            "fraction_pass_stringent_filters": {
                "Description": "Fraction of genotypes which pass stringent hard filters " + stringent_filters,
                "Number": "A",
                "Type": "Float",
            },
            "fraction_pass_medium_filters": {
                "Description": "Fraction of genotypes which pass medium hard filters " + medium_filters,
                "Number": "A",
                "Type": "Float",
            },
            "fraction_pass_relaxed_filters": {
                "Description": "Fraction of genotypes which pass relaxed hard filters " + relaxed_filters,
                "Number": "A",
                "Type": "Float",
            },
            "rf_bin": {
                "Description": "Variant QC random forest bin, model id " + run_hash,
                "Number": "A",
                "Type": "Integer",
            },
            "rf_score": {
                "Description": "Variant QC random forest score, model id " + run_hash,
                "Number": "A",
                "Type": "Float",
            },
            "CSQ": {
                "Description": "Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|UNIPROT_ISOFORM|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|gnomADe_AF|gnomADe_AFR_AF|gnomADe_AMR_AF|gnomADe_ASJ_AF|gnomADe_EAS_AF|gnomADe_FIN_AF|gnomADe_NFE_AF|gnomADe_OTH_AF|gnomADe_SAS_AF|gnomADg_AF|gnomADg_AFR_AF|gnomADg_AMI_AF|gnomADg_AMR_AF|gnomADg_ASJ_AF|gnomADg_EAS_AF|gnomADg_FIN_AF|gnomADg_MID_AF|gnomADg_NFE_AF|gnomADg_OTH_AF|gnomADg_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|PHENOTYPES|Conservation",
                "Number": "A",
                "Type": "String",
            },
            "consequence": {"Description": "Most severe consequence from VEP110", "Number": "A", "Type": "String"},
            "gene": {
                "Description": "Gene affected by the most severe consequence from VEP110",
                "Number": "A",
                "Type": "String",
            },
            "hgnc_id": {
                "Description": "HGNC id of the gene affected by the most severe consequence from VEP110",
                "Number": "A",
                "Type": "String",
            },
        },
        "format": {
            "HetAB": {"Description": "Hetrozygous allele balance", "Number": "A", "Type": "Float"},
            "stringent_filters": {
                "Description": "Pass/fail stringent hard filters " + stringent_filters,
                "Number": "A",
                "Type": "String",
            },
            "medium_filters": {
                "Description": "Pass/fail hard medium filters " + medium_filters,
                "Number": "A",
                "Type": "String",
            },
            "relaxed_filters": {
                "Description": "Pass/fail relaxed hard filters " + relaxed_filters,
                "Number": "A",
                "Type": "String",
            },
        },
    }

    if export_trace:
        # Prepairing step annotation:
        # Extract the analysis_steps from the MatrixTable globals
        analysis_steps = mt.globals.analysis_steps.collect()[0]

        # Prepare the analysis_steps to be included in the VCF header
        analysis_steps_str = [
            f'##analysis_step=<step_name="{step.step_name}", step_description="{step.step_description}">'
            for step in analysis_steps
        ]

        # Write the list to a temporary file to use as part of the VCF header
        filtered_vcf_dir_python = filtered_vcf_dir.replace("file://", "")
        os.makedirs(filtered_vcf_dir_python, exist_ok=True)
        analysis_header_file = os.path.join(filtered_vcf_dir_python, "vcf_analysis_steps_header.txt")
        with open(analysis_header_file, "w") as f:
            for step_str in analysis_steps_str:
                f.write(step_str + "\n")
        trace_header_file_name = "file://" + analysis_header_file
        print(f"=== VCF header written to the file: {analysis_header_file}")
    else:
        trace_header_file_name = None

    # export per chromosome
    chroms = [*[str(i) for i in range(1, 23)], "X", "Y"]
    chromosomes = ["chr" + chr for chr in chroms]
    for chromosome in chromosomes:
        print(f"=== Exporting chromosme {chromosome} to VCF")
        print()
        mt_chrom = mt.annotate_globals(chromosome=chromosome)
        mt_chrom = mt_chrom.filter_rows(mt_chrom.locus.contig == mt_chrom.chromosome)
        outfile = os.path.join(filtered_vcf_dir, f"{chromosome}_hard_filters.vcf.bgz")
        # Export the MatrixTable to a VCF file, including the custom header lines
        hl.export_vcf(mt_chrom, outfile, metadata=metadata, append_to_header=trace_header_file_name)


def main() -> None:
    # set up
    inputs = parse_config()

    data_root: str = inputs["data_root"]
    dataset_name: str = inputs["dataset_name"]

    mtdir = "file://" + os.path.join(data_root, inputs["matrixtables_lustre_dir"])
    runhash = inputs["runhash"]
    hard_filters = inputs["hard_filters"]

    vcf_timestamp = datetime.datetime.now().strftime("%Y%m%d-%H%M")
    filtered_vcf_dir = "file://" + os.path.join(data_root, inputs["vcf_output_dir"], vcf_timestamp)

    # initialise hail
    tmp_dir: str = inputs["tmp_dir"]
    sc = hail_utils.init_hl(tmp_dir)

    mtfile = os.path.join(mtdir, "mt_hard_filter_combinations.mt")
    export_vcfs(mtfile, filtered_vcf_dir, hard_filters, runhash, export_trace=False)
    hail_utils.stop_hl(sc)


if __name__ == "__main__":
    main()
