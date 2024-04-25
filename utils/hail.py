import os
import shutil
import pyspark
import hail as hl


def clear_temp_folder(tmp_dir: str):
    if not tmp_dir.startswith('file://'):
        return
    tmp_dir = tmp_dir.replace('file://', '')
    print(f"=== Cleaning up temporary folder {tmp_dir}")
    for root, dirs, files in os.walk(tmp_dir):
        for f in files:
            os.unlink(os.path.join(root, f))
        for d in dirs:
            shutil.rmtree(os.path.join(root, d))

def init_hl(tmp_dir: str) -> pyspark.SparkContext:
    clear_temp_folder(tmp_dir)
    print("=== Checking for SparkContext ===")
    sc = pyspark.SparkContext()
    # conf = pyspark.SparkConf().setAppName("HailBatch")
    # sc = pyspark.SparkContext.getOrCreate(conf)
    hadoop_config = sc._jsc.hadoopConfiguration()
    print("=== Initializing Hail ===")
    hl.init(sc=sc, tmp_dir=tmp_dir)
    hl.default_reference("GRCh38")
    return sc


def stop_hl(sc: pyspark.SparkContext) -> None:
    hl.stop()
    sc.stop()
