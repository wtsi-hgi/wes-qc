#variant counts per sample for samples sequenced and called at Broad
#after filtering by intervals but no variant QC
#matrixtable was created earlier for Sanger/Broad comparison
import hail as hl
import pyspark
from wes_qc.utils.utils import parse_config
from wes_qc.variant_counts_per_cq_controls.annotation_functions import annotate_cq, annotate_gnomad, get_counts_per_cq, get_median_ca_fraction



def main():
    # set up
    inputs = parse_config()
    mtdir = inputs['matrixtables_lustre_dir']
    resourcedir = inputs['resource_dir']

    # initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    mtfile = mtdir + "mt_after_var_qc_hard_filter_gt.mt"
    mt = hl.read_matrix_table(mtfile)

    #restrict to samples in Sanger
    broad_samples_file = "file:///lustre/scratch123/qc/compare_broad_sanger_vcfs/broad_accs_to_analyse_s.txt"
    broad_sample_ht = hl.import_table(broad_samples_file, delimiter="\t").key_by('s')
    mt = mt.filter_cols(hl.is_defined(broad_sample_ht[mt.s]))

    cqfile = "file:///lustre/scratch123/qc/compare_broad_sanger_vcfs/broad_vcf_samples_in_sanger_filter_by_bait/split_multi_strip_gt/all_consequences.txt"
    
    #cqfile = resourcedir + "all_consequences.txt"
    mtcq = annotate_cq(mt, cqfile)
    gnomad_htfile = resourcedir + "gnomad.exomes.r2.1.1.sites.liftover_grch38.ht"
    mtgnomad = annotate_gnomad(mtcq, gnomad_htfile)

    get_counts_per_cq(mtgnomad)

    get_median_ca_fraction(mtgnomad)
    

if __name__ == '__main__':
    main()
