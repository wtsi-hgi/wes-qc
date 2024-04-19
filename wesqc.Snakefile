# This is the pipeline definition file for WES-QC

import os
import hail as hl

from evaluation import compare_hard_filter_combinations
from utils.hail import init_hl

configfile: "config/inputs.yaml"
configfile: "config/hard-filter-combinations.yaml"

mtdir = config['matrixtables_lustre_dir'].replace("file://","")

rf_dir = config['var_qc_rf_dir'].replace("file://","")
resourcedir = config['resource_dir'].replace("file://","")

plot_dir = config['plots_dir_local'].replace("file://","")
tmp_dir = config["tmp_dir"]

filterCombinationParams = compare_hard_filter_combinations.HardFilterCombinationParams(**config["hard_filters_combinations"])
runhash = filterCombinationParams.runhash
wd = os.path.join(mtdir,runhash)

rule default:
    input:
        outfile_indel = os.path.join(plot_dir,runhash + "_genotype_hard_filter_comparison_indel_5.txt")

rule generate_giab:
    input:
        giab_vcf    = os.path.join(resourcedir, "HG001_GRCh38_benchmark.interval.illumina.vcf.gz"),
        giab_cqfile = os.path.join(resourcedir, "all.interval.illumina.vep.info.txt"),
    output:
        giab_ht_path = directory(os.path.join(mtdir, "giab_annotated.ht"))
    run:
        init_hl(tmp_dir)
        giab_ht = compare_hard_filter_combinations.prepare_giab_ht(
            f"file://{input.giab_vcf}",
            f"file://{input.giab_cqfile}",
            f"file://{mtdir}"
        )
        giab_ht.write(f"file://{output.giab_ht_path}")

rule annotate_mt:
    input:
        mtfile = directory(os.path.join(mtdir,"mt_varqc_splitmulti.mt")),
        rf_htfile = directory(os.path.join(rf_dir,runhash,"_gnomad_score_binning_tmp.ht")),
        cqfile = os.path.join(resourcedir, "all_consequences.txt")
    output:
        mt_annot_path = directory(os.path.join(wd, 'tmp.hard_filters_combs.mt'))
    run:
        init_hl(tmp_dir)
        mt = hl.read_matrix_table(f"file://{input.mtfile}")
        mt = compare_hard_filter_combinations.clean_mt(mt)
        mt_annot = compare_hard_filter_combinations.annotate_with_rf(mt, f"file://{input.rf_htfile}")
        mt_annot = compare_hard_filter_combinations.annotate_cq(mt_annot, f"file://{input.cqfile}")
        mt_annot.write(f"file://{output.mt_annot_path}", overwrite=True)

rule eveluate_filter_combinations:
    input:
        mt_annot_path = rules.annotate_mt.output.mt_annot_path,
        giab_ht_path= rules.generate_giab.output.giab_ht_path,
        pedfile= "/lustre/scratch123/qc/BiB/trios.EGAN.complete.ped",
    output:
        outfile_snv= os.path.join(plot_dir,runhash + "_genotype_hard_filter_comparison_snv_5.txt"),
        outfile_indel=os.path.join(plot_dir,runhash + "_genotype_hard_filter_comparison_indel_5.txt")
    run:
        giab_ht = hl.read_table(f"file://{input.giab_ht_path}")

        results = compare_hard_filter_combinations.filter_and_count(
            mt_path=f"file://{input.mt_annot_path}",
            ht_giab=giab_ht,
            pedfile=f"file://{input.pedfile}",
            mtdir=f"file://{wd}",
            filters=filterCombinationParams
        )

        compare_hard_filter_combinations.print_results(results,output.outfile_snv,'snv')
        compare_hard_filter_combinations.print_results(results,output.outfile_indel,'indel')
