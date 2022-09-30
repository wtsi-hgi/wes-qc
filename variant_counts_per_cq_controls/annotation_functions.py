#functions for annotating with consequence and gnomad and counting variants per cq
import hail as hl
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



def get_counts_per_cq(mt_in: hl.MatrixTable):
    '''
    Get median counts of each consequence per sample
    :param hl.MatrixTable mt_in: Input MatrixTable
    '''
    #split mt by snvs and indels
    snv_mt = mt_in.filter_rows(hl.is_snp(mt_in.alleles[0], mt_in.alleles[1]))
    indel_mt = mt_in.filter_rows(hl.is_indel(mt_in.alleles[0], mt_in.alleles[1]))

    #get median numbers of variants by consequence
    synonymous_counts =  median_count_for_cq(snv_mt, ['synonymous_variant'])
    missense_counts =  median_count_for_cq(snv_mt, ['missense_variant'])
    nonsense_counts =  median_count_for_cq(snv_mt, ['stop_gained'])
    splice_acc_don_counts = median_count_for_cq(snv_mt, ['splice_acceptor_variant', 'splice_donor_variant'])

    frameshift_counts =  median_count_for_cq(indel_mt, ['frameshift_variant'])
    inframe_indel_counts = median_count_for_cq(indel_mt, ['inframe_deletion', 'inframe_insertion'])

    coding_snv_counts = median_count_for_cq(snv_mt, ['synonymous_variant', 'missense_variant', 'stop_gained','splice_acceptor_variant', 'splice_donor_variant', 'sart_lost', 'stop_lost'])
    coding_indel_counts = median_count_for_cq(indel_mt, ['frameshift_variant', 'inframe_deletion', 'inframe_insertion'])

    print("Synonymous: Total " + str(synonymous_counts[0]) + " rare " + str(synonymous_counts[1]))
    print("Missense: Total " + str(missense_counts[0]) + " rare " + str(missense_counts[1]))
    print("Nonsense: Total " + str(nonsense_counts[0]) + " rare " + str(nonsense_counts[1]))
    print("Splicing: Total " + str(splice_acc_don_counts[0]) + " rare " + str(splice_acc_don_counts[1]))
    print("Frameshift: Total " + str(frameshift_counts[0]) + " rare " + str(frameshift_counts[1]))
    print("In-frame indel: Total " + str(inframe_indel_counts[0]) + " rare " + str(inframe_indel_counts[1]))
    print("Coding SNV: Total " + str(coding_snv_counts[0]) + " rare " + str(coding_snv_counts[1]))
    print("Coding indel: Total " + str(coding_indel_counts[0]) + " rare " + str(coding_indel_counts[1]))


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

    print("Median fraction CA per sample total: " + "{:.2f}".format(total_ca_frac))
    print("Median fraction CA per sample rare: " + "{:.2f}".format(rare_ca_frac))


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
