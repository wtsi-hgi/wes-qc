import os.path
import hail as hl
import pyspark
from utils.utils import parse_config


def restrict_to_baits(infile: str, outfile: str, baits: str, padding=50):
    ht = hl.import_locus_intervals(baits)
    ht = ht.annotate(bait=hl.if_else(
        ht.interval.start.contig == 'chrM',
        ht.interval,
        hl.locus_interval(
            contig=ht.interval.start.contig,
            start=ht.interval.start.position-padding,
            end=ht.interval.end.position+padding,
            includes_start=True,
            includes_end=True
        )
    )).key_by('bait')

    mt = hl.read_matrix_table(infile)
    mt_initial = mt.count_rows()
    rows_per_partition = mt_initial // mt.n_partitions()

    mt = mt.filter_rows(hl.is_defined(ht[mt.locus]))
    mt_filtered = mt.count_rows()

    print(f'{mt_initial-mt_filtered} ({(mt_initial-mt_filtered)/mt_initial:.1%}) variants are not in baits')

    mt = mt.repartition(mt_filtered // rows_per_partition)
    mt.write(outfile)


def main():
    inputs = parse_config()
    mtdir = inputs['matrixtables_lustre_dir']

    # initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    mtfile = os.path.join(mtdir, "mt_pops_QC_filters_sample_qc.mt")
    shrinked_file = os.path.join(mtdir, "mt_pops_QC_filters_sample_qc.shrank.mt")
    baits = inputs['baits_file']
    restrict_to_baits(mtfile, shrinked_file, baits)


if __name__ == '__main__':
    main()
