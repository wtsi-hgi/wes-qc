import pyspark
import hail as hl


def init_hl(tmp_dir: str) -> pyspark.SparkContext:
    print("=== Checking for SparkContext ===")
    sc = pyspark.SparkContext()
    #conf = pyspark.SparkConf().setAppName("HailBatch")
    #sc = pyspark.SparkContext.getOrCreate(conf)
    hadoop_config = sc._jsc.hadoopConfiguration()
    print("=== Initializing Hail ===")
    hl.init(sc=sc, tmp_dir=tmp_dir)
    hl.default_reference("GRCh38")
    return sc

def stop_hl(sc: pyspark.SparkContext) -> None:
    hl.stop()
    sc.stop()
