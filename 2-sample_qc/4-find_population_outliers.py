# perform hail sample QC stratified by superpopulation and identify outliers
import hail as hl
import pyspark
from wes_qc.utils.utils import parse_config


def annotate_pops(pruned_mt_file: str, pop_ht_file: str, pops_annot_mt_file: str):
    '''
    Annotate LD pruned matrixtable with predicted populations
    :param str pruned_mt_file: ld pruned MT file
    :param str pop_ht_file: predicted population HT file
    :param str pops_annot_mt_file: population annotated MT file
    '''
    mt = hl.read_matrix_table(pruned_mt_file)
    pop_ht = hl.read_table(pop_ht_file)
    mt = mt.annotate_cols(assigned_pop=pop_ht[mt.s].pop)
    mt.write(pops_annot_mt_file, overwrite=True)



def main():
    #set up
    inputs = parse_config()
    importmtdir = inputs['load_matrixtables_lustre_dir']
    mtdir = inputs['matrixtables_lustre_dir']
    annotdir = inputs['annotation_lustre_dir']

    #initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    #add population labels to ld-pruned biallelic autosomal SNPs
    pruned_mt_file = mtdir + "mt_ldpruned.mt"
    pop_ht_file = mtdir + "pop_assignments.ht"
    pops_annot_mt_file = mtdir + "mt_ldpruned_pops_annotated.mt"
    annotate_pops(pruned_mt_file, pop_ht_file, pops_annot_mt_file)


if __name__ == '__main__':
    main() 