# This is the pipeline definition file for WES-QC

import os
import json
import hail as hl

from evaluation import compare_hard_filter_combinations
from utils.hail import init_hl, stop_hl

configfile: "config/inputs.yaml"
configfile: "config/hard-filter-combinations.yaml"

data_root:str = config['data_root']

mtdir:str = os.path.join(data_root, config['matrixtables_lustre_dir'])

rf_dir:str = os.path.join(data_root, config['var_qc_rf_dir'])
resourcedir:str = os.path.join(data_root, config['resource_dir'])
stats_dir:str = os.path.join(data_root, config['stats_dir'])
plot_dir:str = os.path.join(data_root, config['plots_dir_local'])

tmp_dir:str = config["tmp_dir"]

filterCombinationParams = compare_hard_filter_combinations.HardFilterCombinationParams(**config["hard_filters_combinations"])
runhash:str = filterCombinationParams.runhash
wd = os.path.join(mtdir, runhash)

is_hail_running = False

rule default:
    input:
        outfile_snp = os.path.join(plot_dir, runhash + "_genotype_hard_filter_comparison_snv.txt"),
        outfile_indel = os.path.join(plot_dir, runhash + "_genotype_hard_filter_comparison_indel.txt")

rule generate_giab:
    input:
        giab_vcf    = os.path.join(resourcedir, "HG001_GRCh38_benchmark.interval.illumina.vcf.gz"),
        giab_cqfile = os.path.join(resourcedir, "all.interval.illumina.vep.info.txt"),
    output:
        giab_ht_path = directory(os.path.join(mtdir, "giab_annotated.ht"))
    run:
        sc = init_hl(tmp_dir)
        giab_ht = compare_hard_filter_combinations.prepare_giab_ht(
            f"file://{input.giab_vcf}",
            f"file://{input.giab_cqfile}",
            f"file://{mtdir}"
        )
        giab_ht.write(f"file://{output.giab_ht_path}")
        stop_hl(sc)

rule annotate_mt:
    input:
        mtfile = os.path.join(mtdir,"mt_varqc_splitmulti.mt"),
        rf_htfile = os.path.join(rf_dir,runhash,"_gnomad_score_binning_tmp.ht"),
        cqfile = os.path.join(resourcedir, "all_consequences.txt")
    output:
        mt_annot_path = directory(os.path.join(wd, 'tmp.hard_filters_combs.mt'))
    run:
        sc = init_hl(tmp_dir)
        mt = hl.read_matrix_table(f"file://{input.mtfile}")
        mt = compare_hard_filter_combinations.clean_mt(mt)
        mt_annot = compare_hard_filter_combinations.annotate_with_rf(mt, f"file://{input.rf_htfile}")
        mt_annot = compare_hard_filter_combinations.annotate_cq(mt_annot, f"file://{input.cqfile}")
        mt_annot.write(f"file://{output.mt_annot_path}", overwrite=True)
        stop_hl(sc)

rule filter_var_type:
    input:
        mt_annot_path = rules.annotate_mt.output.mt_annot_path,
    output:
        mt_filtered_path = directory(os.path.join(mtdir,'tmp.hard_filters_combs.{label}.mt'))
    run:
        sc = init_hl(tmp_dir)
        compare_hard_filter_combinations.filter_var_type(
            f"file://{input.mt_annot_path}",
            f"file://{output.mt_filtered_path}",
            wildcards.label
        )
        stop_hl(sc)

rule evaluate_filter_combinations:
    input:
        mt_path = os.path.join(mtdir,"tmp.hard_filters_combs.{label}.mt"),
        giab_ht_path= rules.generate_giab.output.giab_ht_path,
        pedfile= os.path.join(resourcedir, config["pedfile_name"])
    output:
        json_dump_file=os.path.join(wd, 'evaluation.{label}.json')
    run:
        sc = init_hl(tmp_dir)
        giab_ht = hl.read_table(f"file://{input.giab_ht_path}")
        json_dump_folder = os.path.join(stats_dir,'json_dump', runhash)
        os.makedirs(json_dump_folder, exist_ok=True)
        results = compare_hard_filter_combinations.filter_and_count_by_type(
            mt_path=f"file://{input.mt_path}",
            ht_giab=giab_ht,
            pedfile=f"file://{input.pedfile}",
            mtdir=f"file://{wd}",
            filters=filterCombinationParams,
            var_type=wildcards.label,
            json_dump_folder=json_dump_folder,
        )
        with open(output.json_dump_file,'w') as f:
            json.dump(results,f)
        stop_hl(sc)

rule rename_snp_snv: # FIXME: A quich fix to address naming inconsistency
    input:
        json_dump_snp = os.path.join(wd, 'evaluation.snp.json')
    output:
        json_dump_snv = os.path.join(wd, 'evaluation.snv.json')
    shell:
        "ln -s {input} {output}"


rule filter_combination_stats:
    input:
        json_dump_file = os.path.join(wd, 'evaluation.{label}.json')
    output:
        outfile = os.path.join(plot_dir, runhash + "_genotype_hard_filter_comparison_{label}.txt"),
    run:
        with open(input.json_dump_file) as f:
            results = json.load(f)
        compare_hard_filter_combinations.print_results(results,output.outfile,wildcards.label)