"""
The main pipeline definition file for the WES-QC
"""

import os
import json
import hail as hl

from wes_qc import compare_hardfilter
from wes_qc import hail_utils

configfile: "config/inputs.yaml"
configfile: "config/hard-filter-combinations.yaml"

data_root:str = config['data_root']
dataset_name: str = config['dataset_name']
mtdir:str = os.path.join(data_root, config['matrixtables_lustre_dir'])

annot_dir:str = os.path.join(data_root, config['annotation_lustre_dir'])
rf_dir:str = os.path.join(data_root, config['var_qc_rf_dir'])
resourcedir:str = os.path.join(data_root, config['resource_dir'])
stats_dir:str = os.path.join(data_root, config['stats_dir'])
plot_dir:str = os.path.join(data_root, config['plots_dir_local'])

tmp_dir:str = config["tmp_dir"]

filterCombinationParams = compare_hardfilter.HardFilterCombinationParams(**config["hard_filters_combinations"])
runhash:str = filterCombinationParams.runhash
wd = os.path.join(mtdir, runhash)

is_hail_running = False

rule default:
    input:
        out_file_snv=os.path.join(wd, 'evaluation.snv.json'),
        #out_file_indel=os.path.join(wd, 'evaluation.indel.json')
        #outfile_snp = os.path.join(plot_dir, runhash, "evaluation.snv.csv"),
        #outfile_indel = os.path.join(plot_dir, runhash, "evaluation.indel.csv")

rule generate_giab:
    input:
        giab_vcf    = os.path.join(resourcedir, "HG001_GRCh38_benchmark.interval.illumina.vcf.gz"),
        giab_cqfile = os.path.join(resourcedir, "all.interval.illumina.vep.info.txt"),
    output:
        giab_ht_path = directory(os.path.join(mtdir, "giab_annotated.ht"))
    benchmark: os.path.join(stats_dir, 'benchmark_generate_giab.txt')
    run:
        sc = hail_utils.init_hl(tmp_dir)
        giab_ht = compare_hardfilter.prepare_giab_ht(
            f"file://{input.giab_vcf}",
            f"file://{input.giab_cqfile}",
            f"file://{mtdir}"
        )
        giab_ht.write(f"file://{output.giab_ht_path}")
        hail_utils.stop_hl(sc)

rule annotate_mt:
    input:
        mtfile = os.path.join(mtdir,"mt_varqc_splitmulti.mt"),
        rf_htfile = os.path.join(rf_dir,runhash,"_gnomad_score_binning_tmp.ht"),
        cqfile = os.path.join(resourcedir, "all_consequences.txt")
    output:
        mt_annot_path = directory(os.path.join(wd, 'tmp.hard_filters_combs.mt'))
    benchmark: os.path.join(stats_dir,'benchmark_annotate_mt.txt')
    run:
        sc = hail_utils.init_hl(tmp_dir)
        mt = hl.read_matrix_table(f"file://{input.mtfile}")
        mt = compare_hardfilter.clean_mt(mt)
        mt_annot = compare_hardfilter.annotate_with_rf(mt,f"file://{input.rf_htfile}")
        mt_annot = compare_hardfilter.annotate_cq(mt_annot,f"file://{input.cqfile}")
        mt_annot.write(f"file://{output.mt_annot_path}", overwrite=True)
        hail_utils.stop_hl(sc)

rule filter_var_type:
    input:
        mt_annot_path = rules.annotate_mt.output.mt_annot_path,
    output:
        mt_filtered_path = directory(os.path.join(mtdir,'tmp.hard_filters_combs.{label}.mt'))
    benchmark: os.path.join(stats_dir,'benchmark_{label}_filter_var_type.txt')
    run:
        sc = hail_utils.init_hl(tmp_dir)
        compare_hardfilter.filter_var_type(
            f"file://{input.mt_annot_path}",
            f"file://{output.mt_filtered_path}",
            wildcards.label
        )
        hail_utils.stop_hl(sc)

rule evaluate_filter_combinations:
    input:
        mt_path = os.path.join(mtdir,"tmp.hard_filters_combs.{label}.mt"),
        giab_ht_path= rules.generate_giab.output.giab_ht_path,
        pedfile= os.path.join(annot_dir, config["pedfile_name"])
    output:
        json_dump_file=protected(os.path.join(wd, 'evaluation.{label}.json'))
    benchmark: os.path.join(stats_dir,'benchmark_{label}_evaluate_filter_combinations.txt')
    run:
        sc = hail_utils.init_hl(tmp_dir)
        giab_ht = hl.read_table(f"file://{input.giab_ht_path}")
        json_dump_folder = os.path.join(stats_dir,'json_dump', runhash)
        os.makedirs(json_dump_folder, exist_ok=True)
        results = compare_hardfilter.filter_and_count_by_type(
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
        hail_utils.stop_hl(sc)

rule rename_snp_snv: # FIXME: A quick fix to address naming inconsistency
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
        #outfile = os.path.join(wd, 'evaluation.{label}.csv'),
        outfile = os.path.join(plot_dir, runhash, "evaluation.{label}.csv")
    # For the BiB version we use the R script to make final graphs and tables instead of Python code
    params:
        var_type = "{label}",
        basename = 'evaluation.{label}',
        target_dir = os.path.join(plot_dir, runhash)
    shell:
        """Rscript scripts/plot_hardfilter.r {input} {params.var_type}
        mkdir -p {params.target_dir}
        cp {params.basename}.csv {params.target_dir}
        cp {params.basename}.*.png {params.target_dir}
        """

#    run:
#        with open(input.json_dump_file) as f:
#            results = json.load(f)
#        compare_hardfilter.print_results(results,output.outfile,wildcards.label)
