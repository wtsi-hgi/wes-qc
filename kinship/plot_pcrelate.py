# create plots from binned random forest output
import hail as hl
import pyspark
import os
from bokeh.plotting import output_file, save
from wes_qc.utils.utils import parse_config


def plot_kinship(pcrelate_htfile: str, pedfile: str, annotfile: str, plot_dir: str):
    '''
    :param str pcrelate_htfile: hail table output file from pc-relate
    :param str pedfile: pedfile
    :param str plot_dir: Direcotry for output plot
    '''
    ht = hl.read_table(pcrelate_htfile)
    ht = ht.annotate(pair = (ht.i.s + "_" + ht.j.s))

    fam = hl.import_fam(pedfile)
    probands = fam.id.collect()
    mums = fam.mat_id.collect()
    dads = fam.pat_id.collect()
    fam_members = probands + mums + dads

    #only want to keep the samples in trios
    ht = ht.filter( (hl.set(fam_members).contains(ht.i.s)) & (hl.set(fam_members).contains(ht.j.s)), keep = True)
    print("Filtered by sample:")
    ht.count()

    #annnotate with pair id
    ht = ht.annotate(pair = (ht.i.s + "_" + ht.j.s))

    #create a list of all possible parent_child pairs (both orientations)
    parent_child1 = [m + "_" + n for m,n in zip(probands,mums)]
    parent_child2 = [m + "_" + n for m,n in zip(probands,dads)]
    parent_child3 = [m + "_" + n for m,n in zip(mums,probands)]
    parent_child4 = [m + "_" + n for m,n in zip(dads,probands)]
    parent_child_all = parent_child1 + parent_child2 + parent_child3 + parent_child4

    #annotate with family relationships
    annotexpr = hl.case().when(hl.set(parent_child_all).contains(ht.pair), "parent_child").default("unrelated")
    ht = ht.annotate(relationship = annotexpr)

    ht.write(annotfile, overwrite=True)

    #plot
    plotfile = plot_dir + "pc-relate.html"
    p = hl.plot.scatter(ht.ibd0, ht.kin, label = ht.relationship, xlabel='IBD0', ylabel='kin', title = "PC-relate for samples in trios")
    output_file(filename=plotfile)
    save(p)


def main():
    #set up
    inputs = parse_config()
    mtdir = inputs['matrixtables_lustre_dir']
    resourcedir = inputs['resource_dir']
    root_plot_dir = inputs['plots_dir_local']

    #initialise hail
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    pcrelate_htfile = mtdir + "pcrelate.ht"
    annotfile = mtdir +  "pcrelate_annotated.ht"
    pedfile = resourcedir + "trios.ped"
    plot_dir = root_plot_dir + "kinship/"

    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    plot_kinship(pcrelate_htfile, pedfile, annotfile, plot_dir)

if __name__ == '__main__':
    main() 