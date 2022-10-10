#functions for annotating with consequence and gnomad and counting variants per cq
import hail as hl
import os
import numpy as np

def annotate_cq(mt: hl.MatrixTable, cqfile: str) -> hl.MatrixTable:
    '''
    Annotate a hail MatrixTable from a tsv file of consequence annotations
    :param hl.MatrixTable mt: input MatrixTable
    :param str cqfile: path to tsv file containing consequence annotation
    '''
    ht=hl.import_table(cqfile,types={'f0':'str','f1':'int32', 'f2':'str','f3':'str','f4':'str', 'f5':'str'}, no_header=True)
    ht=ht.annotate(chr=ht.f0)
    ht=ht.annotate(pos=ht.f1)
    ht=ht.annotate(rs=ht.f2)
    ht=ht.annotate(ref=ht.f3)
    ht=ht.annotate(alt=ht.f4)
    ht=ht.annotate(consequence=ht.f5)
    ht = ht.key_by(
    locus=hl.locus(ht.chr, ht.pos), alleles=[ht.ref,ht.alt])
    ht=ht.drop(ht.f0,ht.f1,ht.f2,ht.f3,ht.f4,ht.chr,ht.pos,ht.ref,ht.alt)
    ht = ht.key_by(ht.locus, ht.alleles)

    mt=mt.annotate_rows(consequence=ht[mt.row_key].consequence)
    return mt


def annotate_gnomad(mt_in: hl.MatrixTable, gnomad_htfile: str) -> hl.MatrixTable:
    '''
    Annotate matrixtable with AC and AF from gnomad
    :param hl.MatrixTable mt_in: Input MatrixTable
    :param str gnomad_htfile: gnomAD Hail table file
    :return: Annotated MatrixTable
    '''
    gnomad_ht = hl.read_table(gnomad_htfile)
    #ac = gnomadht.freq[0].AC
    #af = gnomadht.freq[0].AF
    mt_in = mt_in.annotate_rows(gnomad_AF=gnomad_ht[mt_in.row_key].freq[0].AF)
    mt_in = mt_in.annotate_rows(gnomad_AC=gnomad_ht[mt_in.row_key].freq[0].AC)

    return mt_in



def get_counts_per_cq(mt_in: hl.MatrixTable, outfile: str):
    '''
    Get median counts of each consequence per sample
    :param hl.MatrixTable mt_in: Input MatrixTable
    :param str outfile: Output file path
    '''
    #split mt by snvs and indels
    snv_mt = mt_in.filter_rows(hl.is_snp(mt_in.alleles[0], mt_in.alleles[1]))
    indel_mt = mt_in.filter_rows(hl.is_indel(mt_in.alleles[0], mt_in.alleles[1]))

    #get median numbers of variants by consequence
    # synonymous_counts =  median_count_for_cq(snv_mt, ['synonymous_variant'])
    # missense_counts =  median_count_for_cq(snv_mt, ['missense_variant'])
    # nonsense_counts =  median_count_for_cq(snv_mt, ['stop_gained'])
    # splice_acc_don_counts = median_count_for_cq(snv_mt, ['splice_acceptor_variant', 'splice_donor_variant'])

    # frameshift_counts =  median_count_for_cq(indel_mt, ['frameshift_variant'])
    # inframe_indel_counts = median_count_for_cq(indel_mt, ['inframe_deletion', 'inframe_insertion'])

    # coding_snv_counts = median_count_for_cq(snv_mt, ['synonymous_variant', 'missense_variant', 'stop_gained','splice_acceptor_variant', 'splice_donor_variant', 'sart_lost', 'stop_lost'])
    # coding_indel_counts = median_count_for_cq(indel_mt, ['frameshift_variant', 'inframe_deletion', 'inframe_insertion'])

    # print("Synonymous: Total " + str(synonymous_counts[0]) + " rare " + str(synonymous_counts[1]))
    # print("Missense: Total " + str(missense_counts[0]) + " rare " + str(missense_counts[1]))
    # print("Nonsense: Total " + str(nonsense_counts[0]) + " rare " + str(nonsense_counts[1]))
    # print("Splicing: Total " + str(splice_acc_don_counts[0]) + " rare " + str(splice_acc_don_counts[1]))
    # print("Frameshift: Total " + str(frameshift_counts[0]) + " rare " + str(frameshift_counts[1]))
    # print("In-frame indel: Total " + str(inframe_indel_counts[0]) + " rare " + str(inframe_indel_counts[1]))
    # print("Coding SNV: Total " + str(coding_snv_counts[0]) + " rare " + str(coding_snv_counts[1]))
    # print("Coding indel: Total " + str(coding_indel_counts[0]) + " rare " + str(coding_indel_counts[1]))

    # get variant counts per consequence
    synonymous_all, synonymous_rare = counts_per_cq(snv_mt, ['synonymous_variant'])
    missense_all, missense_rare = counts_per_cq(snv_mt, ['missense_variant'])
    nonsense_all, nonsense_rare = counts_per_cq(snv_mt, ['stop_gained'])
    splice_all, splice_rare = counts_per_cq(snv_mt, ['splice_acceptor_variant', 'splice_donor_variant'])

    frameshift_all, frameshift_rare = counts_per_cq(indel_mt, ['frameshift_variant'])
    inframe_all, inframe_rare = counts_per_cq(indel_mt, ['inframe_deletion', 'inframe_insertion'])

    coding_snv_all, coding_snv_rare = counts_per_cq(
        snv_mt, ['synonymous_variant', 'missense_variant', 'stop_gained', 'splice_acceptor_variant', 
        'splice_donor_variant', 'sart_lost', 'stop_lost'])
    coding_indel_all, coding_indel_rare = counts_per_cq(
        indel_mt, ['frameshift_variant', 'inframe_deletion', 'inframe_insertion'])

    # print out medians
    print ("Median variant counts per sample for each functional class:")
    print("Synonymous: All variants" + str(np.median(synonymous_all)) + " rare " + str(np.median(synonymous_rare)))
    print("Missense: All variants " + str(np.median(missense_all)) + " rare " + str(np.median(missense_rare)))
    print("Nonsense: All variants " + str(np.median(nonsense_all)) + " rare " + str(np.median(nonsense_rare)))
    print("Splicing: All variants " + str(np.median(splice_all)) + " rare " + str(np.median(splice_rare)))
    print("Frameshift: All variants " + str(np.median(frameshift_all)) + " rare " + str(np.median(frameshift_rare)))
    print("In-frame indel: All variants " + str(np.median(inframe_all)) + " rare " + str(np.median(inframe_rare)))
    print("Coding SNV: All variants " + str(np.median(coding_snv_all)) + " rare " + str(np.median(coding_snv_rare)))
    print("Coding indel: All variants " + str(np.median(coding_indel_all)) + " rare " + str(np.median(coding_indel_rare)))


    # print to output table
    outdata = list(zip(synonymous_all, synonymous_rare, missense_all, missense_rare, nonsense_all, nonsense_rare, 
        splice_all, splice_rare, frameshift_all, frameshift_rare, inframe_all, inframe_rare, coding_snv_all, 
        coding_snv_rare, coding_indel_all, coding_indel_rare))

    header = ['synonymous_all', 'synonymous_rare', 'missense_all', 'missense_rare', 'nonsense_all', 
        'nonsense_rare', 'splice_all', 'splice_rare', 'frameshift_all', 'frameshife_rare', 'in_frame_all', 
        'in_frame_rare', 'coding_snv_all', 'coding_snv_rare', 'coding_indel_all', 'coding_indel_rare']

    n_rows = len(synonymous_all)
    n_cols = len(header)
    with open(outfile, 'w') as o:
        o.write(("\t").join(header))
        o.write("\n")
        for i in range(0, n_rows):
            linedata = []
            for j in range(0, n_cols):
                linedata.append(str(outdata[i][j]))
            o.write(("\t").join(linedata))
            o.write("\n")


def counts_per_cq(mt_in: hl.MatrixTable, cqs: list) -> tuple:
    '''
    Get all counts per sample for a list of consequences
    :param hl.MatrixTable mt_in: Input MatrixTable
    :param list cqs: List fof consequences
    :return: dict
    '''
    mt = mt_in.filter_rows(hl.literal(cqs).contains(mt_in.info.consequence))
    mt_rare = mt.filter_rows(mt.gnomad_AC < 5)
    mt = hl.sample_qc(mt)
    mt_rare = hl.sample_qc(mt_rare)
    sampleqc_ht = mt.cols()
    sampleqc_rare_ht = mt_rare.cols()
    counts_all = sampleqc_ht.sample_qc.n_non_ref.collect()
    counts_rare = sampleqc_rare_ht.sample_qc.n_non_ref.collect()

    #remove aggregate intermediates from tmp - this is a hack as this dir fills up and causes this to exit
    aggdir = "/lustre/scratch123/qc/tmp/aggregate_intermediates/"
    aggfiles = os.listdir(aggdir)
    for af in aggfiles:
        aggpath = aggdir + af
        os.remove(aggpath)

    return counts_all, counts_rare


def median_count_for_cq(mt_in: hl.MatrixTable, cqs: list) -> tuple:
    '''
    Get median counts per sample for a list of consequences
    :param hl.MatrixTable mt_in: Input MatrixTable
    :param list cqs: List fof consequences
    :return: tuple
    '''

    mt = mt_in.filter_rows(hl.literal(cqs).contains(mt_in.consequence))
    mt_rare = mt.filter_rows(mt.gnomad_AC < 5)

    mt = hl.sample_qc(mt)
    mt_rare = hl.sample_qc(mt_rare)
    sampleqc_ht = mt.cols()
    sampleqc_rare_ht = mt_rare.cols()
    total_median = hl.median(sampleqc_ht.sample_qc.n_non_ref.collect()).collect()[0]
    rare_median = hl.median(sampleqc_rare_ht.sample_qc.n_non_ref.collect()).collect()[0]
    
    return total_median, rare_median


def get_ca_per_sample(mt_in: hl.MatrixTable) -> list:
    '''
    Get fraction CA for each sample in a matrixtable
    :param hl.MatrixTable mt_in: Input MatrixTable
    :return: list of C>A fractions
    '''
    mt_in = mt_in.annotate_rows(is_CA=((mt_in.alleles[0] == "C") & (mt_in.alleles[1] == "A")) | (
        (mt_in.alleles[0] == "G") & (mt_in.alleles[1] == "T")))
    ca_mt = mt_in.filter_rows(mt_in.is_CA == True)
    # run sample qc and extract cols
    mt_in = hl.sample_qc(mt_in)
    ca_mt = hl.sample_qc(ca_mt)

    #remove aggregate intermediates from tmp - this is a hack as this dir fills up and causes this to exit
    aggdir = "/lustre/scratch123/qc/tmp/aggregate_intermediates/"
    aggfiles = os.listdir(aggdir)
    for af in aggfiles:
        aggpath = aggdir + af
        os.remove(aggpath)

    snv_qc = mt_in.cols()
    ca_qc = ca_mt.cols()

    snv_count = snv_qc.sample_qc.n_non_ref.collect()
    ca_count = ca_qc.sample_qc.n_non_ref.collect()
    fracs = []
    for i in range (0, len(snv_count)):
        if snv_count[i] > 0:
            f = ca_count[i]/snv_count[i]
            fracs.append(f)
        else:
            fracs.append(0.0)

    return fracs


def get_ca_fractions(mt_in: hl.MatrixTable, outfile: str):
    '''
    Get fraction CA per sample
    :param hl.MatrixTable mt_in: Input MatrixTable
    :param str outfile: output text file path
    '''
    # filter to SNVs and to rare
    snv_mt = mt_in.filter_rows(hl.is_snp(mt_in.alleles[0], mt_in.alleles[1]))
    snv_mt_rare = snv_mt.filter_rows(snv_mt.gnomad_AC < 5)

    # total_ca_frac = get_median_ca_per_sample(snv_mt)
    # rare_ca_frac = get_median_ca_per_sample(snv_mt_rare)
    total_ca_frac = get_ca_per_sample(snv_mt)
    rare_ca_frac = get_ca_per_sample(snv_mt_rare)

    print("Median fraction CA per sample total: " + "{:.5f}".format(np.median(total_ca_frac)))
    print("Median fraction CA per sample rare: " + "{:.5f}".format(np.median(rare_ca_frac)))
    outdata = list(zip(total_ca_frac, rare_ca_frac))

    header = ['total_ca_fraction', 'rare_ca_fraction']
    n_rows = len(outdata)
    n_cols = len(header)

    with open(outfile, 'w') as o:
        o.write(("\t").join(header))
        o.write("\n")
        for i in range(0, n_rows):
            linedata = []
            for j in range(0, n_cols):
                n = "{:.5f}".format(outdata[i][j])
                linedata.append(n)
            o.write(("\t").join(linedata))
            o.write("\n")


def get_median_ca_fraction(mt_in: hl.MatrixTable):
    '''
    Get fraction CA per sample
    :param hl.MatrixTable mt_in: Input MatrixTable
    '''
    #filter to SNVs and to rare
    snv_mt = mt_in.filter_rows(hl.is_snp(mt_in.alleles[0], mt_in.alleles[1]))
    snv_mt_rare = snv_mt.filter_rows(snv_mt.gnomad_AC < 5)

    total_ca_frac = get_median_ca_per_sample(snv_mt)
    rare_ca_frac = get_median_ca_per_sample(snv_mt_rare)

    print("Median fraction CA per sample total: " + "{:.5f}".format(total_ca_frac))
    print("Median fraction CA per sample rare: " + "{:.5f}".format(rare_ca_frac))


def get_median_ca_per_sample(mt_in: hl.MatrixTable) -> float:
    '''
    Get fraction CA per sample
    :param hl.MatrixTable mt_in: Input MatrixTable
    :return: float
    ''' 
    #annotate SNVs with is_CA
    mt_in = mt_in.annotate_rows(is_CA=((mt_in.alleles[0] == "C") & (mt_in.alleles[1] == "A")) | ((mt_in.alleles[0] == "G") & (mt_in.alleles[1] == "T")))
    ca_mt = mt_in.filter_rows(mt_in.is_CA==True)
    #run sample qc and extract cols
    mt_in = hl.sample_qc(mt_in)
    ca_mt = hl.sample_qc(ca_mt)
    snv_qc = mt_in.cols()
    ca_qc = ca_mt.cols()

    snv_samples = snv_qc.s.collect()
    ca_samples = ca_qc.s.collect()
    snv_count = snv_qc.sample_qc.n_non_ref.collect()
    ca_count = ca_qc.sample_qc.n_non_ref.collect()
    snv_count_sample = {snv_samples[i]: snv_count[i] for i in range(len(snv_samples))}
    ca_count_sample = {ca_samples[i]: ca_count[i] for i in range(len(ca_samples))}

    sample_fracs = []
    for s in snv_count_sample.keys():
        tot_vars = int(snv_count_sample[s])
        ca_vars = int(ca_count_sample[s])
        if tot_vars > 0:
            s_frac = ca_vars/tot_vars
            sample_fracs.append(s_frac)
        else:
            sample_fracs.append(0.0)

    median_ca = np.median(sample_fracs)

    return median_ca
