#create a matrixtable containing only broad samples called at Sanger for variant counts
import hail as hl
import pyspark
from wes_qc.utils.utils import parse_config


def main():
    # set up
    inputs = parse_config()
    mtdir = inputs['matrixtables_lustre_dir']

    # initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    input_mtfile = "file:///lustre/scratch123/qc/matrixtables/gatk_unprocessed_with_pop_and_runid.mt"
    output_mtfile = mtdir + "gatk_broad_crams_sanger_calls.mt"
    mt = hl.read_matrix_table(input_mtfile)
    mt = mt.filter_cols(mt.sequencing_location == 'Bristol')#crams came from Bristol but were sequenced at Broad
    mt.write(output_mtfile, overwrite=True)


if __name__ == '__main__':
    main()
