#run transmission disequilibrium test to calculate number of transmitted/untransmitted variants
#annotate with gnomad MAF

import hail as hl
import pyspark
from wes_qc.utils.utils import parse_config


def run_tdt(varqc_mtfile, pedfile, gnomad_htfile, tdt_htfile):
    '''
    Run transmission disequilibrium test, annotate and save output
    '''
    mt = hl.read_matrix_table(varqc_mtfile)
    pedigree = hl.Pedigree.read(pedfile)
    tdt_ht = hl.transmission_disequilibrium_test(mt, pedigree)
    #add gnomad AF
    gnomad_ht = hl.read_table(gnomad_htfile)
    tdt_ht=tdt_ht.annotate(gnomad_af=gnomad_ht[tdt_ht.key].maf)
    #save ht to file
    tdt_ht.write(tdt_htfile, overwrite=True)


def main():
    #set up input variables
    inputs = parse_config()
    mtdir = inputs['matrixtables_lustre_dir']
    resourcedir = inputs['resource_dir']

    #initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    varqc_mtfile = mtdir + "/mts_from_40eedd2d/mt_varqc_splitmulti.mt"#use file from a run of variant QC where no variants exlcuded
    pedfile = resourcedir + "trios.ped"
    gnomad_htfile = resourcedir + "gnomad_v3-0_AF.ht"
    tdt_htfile = mtdir + "tansmission_disequilibrium_test.ht"
    run_tdt(varqc_mtfile, pedfile, gnomad_htfile, tdt_htfile)


if __name__ == '__main__':
    main() 