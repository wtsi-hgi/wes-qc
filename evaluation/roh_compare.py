import hail as hl
import pyspark
import pandas
from utils.utils import rm_mt


def read_roh(path: str):
    roh = pandas.read_csv(path, sep='\s+',
                          usecols=['FID', 'CHR', 'POS1', 'POS2', 'NSNP', 'PHET'],
                          dtype={'FID': str, 'CHR': str})
    roh['NHET'] = (roh['NSNP'] * roh['PHET']).round().astype(int)
    roh['POS_LEN'] = roh['POS2'] - roh['POS1']
    return roh


def filter_roh(df: pandas.DataFrame, min_length_bp=5e6, trim_bp=2e6) -> pandas.DataFrame:
    # ROH length filter
    count_init = df.shape[0]

    df = df[df['POS_LEN'] > min_length_bp]
    count_2 = df.shape[0]

    # exclude regions intersecting with HLA
    hla = ('6', 29941260, 29949572)
    df = df[~((df['CHR'] == hla[0]) &
              (df['POS1'] < hla[2]) &
              (df['POS2'] > hla[1]))]
    count_3 = df.shape[0]

    # remove ROHs with hets
    df = df[df['NHET'] < 1]
    count_4 = df.shape[0]

    print(f'There are {count_init} ROHs\n'
          f'{count_init-count_2} ROHs are less {min_length_bp/1e6} Mb\n'
          f'{count_2-count_3} intersect with HLA\n'
          f'{count_3-count_4} have >1 hets')

    # trim ROHs from each side
    df['POS1'] = df['POS1'] + int(trim_bp)
    df['POS2'] = df['POS2'] - int(trim_bp)
    df['POS_LEN'] = df['POS_LEN'] - 2*int(trim_bp)
    df.drop(['NSNP'], axis=1, inplace=True)

    return df


def prepare_roh(df: pandas.DataFrame, id_map_path: str = None) -> pandas.DataFrame:
    id_map = pandas.read_csv(id_map_path, sep='\t', dtype=str)

    df = df.merge(id_map, left_on='FID', right_on='OrageneID.GSA')
    df.drop(['OrageneID.GSA', 'FID'], axis=1, inplace=True)
    df.rename(columns={'OrageneID.WES': 'ID'}, inplace=True)
    df['CHR'] = 'chr' + df['CHR']

    return df


def make_roh(path: str, map_path: str) -> str:

    roh = read_roh(path)
    roh = filter_roh(roh)
    roh = prepare_roh(roh, id_map_path=map_path)

    new_path = path.replace('.hom', '.hail-ready.hom')
    roh.to_csv(new_path, sep='\t', index=False)

    return new_path


def import_roh(path: str) -> hl.Table:
    roh = hl.import_table(path, key='ID', impute=True, types={'CHR': hl.tstr})
    roh = roh.annotate(region=hl.locus_interval(roh.CHR, roh.POS1, roh.POS2))
    # roh = roh.annotate(locus=hl.struct(chrom=roh.CHR, start=roh.POS1, end=roh.POS2))
    roh = roh.select(roh.region)
    return roh.cache()


def count_hets_in_rohs(mt_path: str, roh_path: str, temp_path: str):

    mt = hl.read_matrix_table(mt_path)
    roh = import_roh(roh_path)
    # samples = roh.aggregate(hl.agg.collect_as_set(roh.ID))

    # print(f'there are {mt.count()} variants in orig mt')
    # data = mt.annotate_rows(roh=roh[mt.locus]).entries()
    # print(f'\nthere are {data.count()} variants after join')
    # data = data.filter(data.s == data.roh.ID)
    # print(f'\nthere are {data.count()} variants after filter')
    # data = data.group_by(data.s, data.roh.locus).aggregate(nhet=hl.agg.count_where(data['GT'].is_het()))

    gt = mt.entries().key_by('s')
    data = gt.join(roh, how='inner')

    data = data.filter(data.region.contains(data.locus))
    data = data.group_by(data.s, data.region).aggregate(nhet=hl.agg.count_where(data['GT'].is_het()),
                                                        nsnp=hl.agg.count())

    data.write(temp_path, overwrite=True)

    # for sample in samples:
    #     sample_roh = roh.filter(roh.ID == sample).key_by('region')
    #
    #     sample_mt = mt.filter_cols(mt.s == sample)
    #     sample_mt = sample_mt.annotate_rows(roh=sample_roh[sample_mt.locus].locus)
    #     sample_mt = sample_mt.filter_rows(hl.is_defined(sample_mt.roh))
    #     # sample_mt = sample_mt.checkpoint(temp_path, overwrite=True)
    #
    #     data = sample_mt.group_rows_by(sample_mt.roh).aggregate(nhet=hl.agg.count_where(sample_mt['GT'].is_het()))
    #     data.nhet.show()
    #
    #     # from gnomad.sample_qc.filtering import compute_stratified_sample_qc
    #     # data = compute_stratified_sample_qc(sample_mt, strata={'roh': sample_mt.roh}, tmp_ht_prefix=None)
    #     #
    #     # data = hl.variant_qc(sample_mt, name='sample_variant_stats')
    #     # sample_mt = sample_mt.annotate_rows(nhet=data.index_rows(sample_mt.row_key).sample_variant_stats.n_het)
    #     # data.group_rows_by(data.roh).aggregate(sum_het=hl.agg.sum(data.sample_variant_stats.n_het))
    #     pass
    #
    # rm_mt(temp_path)


def shrink_input_mt(mt_path: str, roh_path: str):
    mt = hl.read_matrix_table(mt_path)
    roh = import_roh(roh_path).key_by('region')
    samples = roh.aggregate(hl.agg.collect_as_set(roh.ID))

    mt = mt.select_entries(mt.GT)
    mt = mt.drop(mt.assigned_pop, *mt.row_value)
    mt = mt.filter_cols(hl.set(samples).contains(mt.s))
    mt = mt.filter_rows(hl.is_defined(roh[mt.locus]))

    new_path = mt_path.replace('.mt', '.filtered.mt')
    mt.write(new_path)
    return new_path


def main():
    tmp_dir = "hdfs://spark-master:9820/"
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    path = '/lustre/scratch123/projects/gnh_industry/Genes_and_Health_2023_02_44k/GH_44k_autosome_maf0.01_geno0.01_hwe1e-6_ROH_CALLING_OUT.hom'
    map_path = '/lustre/scratch123/projects/gnh_industry/Genes_and_Health_2023_02_44k/link_OrageneID_all-WES_GSA.txt'
    mt_path = 'file:///lustre/scratch123/projects/gnh_industry/Genes_and_Health_2023_02_44k/qc/matrixtables/mt_after_var_qc.50809bee-77-49.mt'
    temp_path = 'file:///lustre/scratch123/projects/gnh_industry/Genes_and_Health_2023_02_44k/qc/matrixtables/mt_after_var_qc.50809bee-77-49.roh-stat.ht'

    roh_path = '/lustre/scratch123/projects/gnh_industry/Genes_and_Health_2023_02_44k/GH_44k_autosome_maf0.01_geno0.01_hwe1e-6_ROH_CALLING_OUT.hail-ready.hom'
    mt_shrinked_path = 'file:///lustre/scratch123/projects/gnh_industry/Genes_and_Health_2023_02_44k/qc/matrixtables/mt_after_var_qc.50809bee-77-49.filtered.mt'
    # roh_path = make_roh(path=path, map_path=map_path)
    # mt_shrinked_path = shrink_input_mt(mt_path, roh_path='file://'+roh_path)

    count_hets_in_rohs(mt_path=mt_shrinked_path, roh_path='file://'+roh_path, temp_path=temp_path)


if __name__ == '__main__':
    main()
