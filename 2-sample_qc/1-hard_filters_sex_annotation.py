#apply gnomad's hard filters and impute sex
#input gatk_unprocessed.mt from step 1.1
import hail as hl
import pyspark
from wes_qc.utils.utils import parse_config

def apply_hard_filters(mt: hl.MatrixTable, mtdir: str) -> hl.MatrixTable:
    '''
    Applies hard filters and annotates samples in the filtered set with call rate
    :param MatrixTable mt: MT containing samples to be ascertained for sex
    :param str mtdir: directory output matrix tables are written to
    :return: MatrixTable with hard filtering annotation
    :rtype: MatrixTable
    '''
    print("Applying hard filters")
    filtered_mt_file = mtdir + "mt_hard_filters_annotated.mt"
    mt = mt.filter_rows((hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1]) &
        (hl.agg.mean(mt.GT.n_alt_alleles()) / 2 > 0.001) &
        (hl.agg.fraction(hl.is_defined(mt.GT)) > 0.99))
    mt = mt.annotate_cols(callrate=hl.agg.fraction(hl.is_defined(mt.GT)))
    mt.write(filtered_mt_file, overwrite=True)

    return mt


def impute_sex(mt: hl.MatrixTable, mtdir: str, annotdir: str, male_threshold: float = 0.8, female_threshold: float = 0.5) -> hl.MatrixTable:
    '''
    Imputes sex, exports data, and annotates mt with this data
    :param MatrixTable mt: MT containing samples to be ascertained for sex
    :param str mtdir: directory output matrix tables are written to
    :param str annotdir: directory annotation files are written to
    :return: MatrixTable with imputed sex annotations stashed in column annotation 'sex_check'
    :rtype: MatrixTable
    '''
    print("Imputing sex with male_threshold = " + str(male_threshold) + " and female threshold = " + str(female_threshold))

    #filter to X and select unphased diploid genotypes
    mt1 = hl.filter_intervals(mt, [hl.parse_locus_interval('chrX')])
    mtx_unphased = mt1.select_entries(GT=hl.unphased_diploid_gt_index_call(mt1.GT.n_alt_alleles()))
    #imput sex on the unphased diploid GTs
    sex_ht = hl.impute_sex(mtx_unphased.GT, aaf_threshold=0.05, female_threshold=female_threshold, male_threshold=male_threshold)
    #export
    sex_ht.export(annotdir + '/sex_annotated.sex_check.txt.bgz')
    #annotate input (all chroms) mt with imputed sex and write to file
    sex_colnames = ['f_stat', 'is_female']
    sex_ht = sex_ht.select(*sex_colnames)
    mt = mt.annotate_cols(**sex_ht[mt.col_key])
    sex_mt_file = mtdir + "mt_sex_annotated.mt"
    print("Writing to " + sex_mt_file)
    mt.write(sex_mt_file, overwrite=True)

    return mt


def annotate_ambiguous_sex(mt: hl.MatrixTable, mtdir: str):
    '''
    Annotated samples with ambiguous sex after impute_sex has run and save as a Table
    :param MatrixTable mt: MT containing imputed sex in column 'sex_check'
    :param str mtdir: directory output matrix tables are written to
    '''
    print("Annotating samples with ambiguous sex:")
    qc_ht = mt.cols()

    qc_ht = qc_ht.annotate(ambiguous_sex=((qc_ht.f_stat >= 0.5) & (hl.is_defined(qc_ht.normalized_y_coverage) & (qc_ht.normalized_y_coverage <= 0.1))) |
                           (hl.is_missing(qc_ht.f_stat)) |
                           ((qc_ht.f_stat >= 0.4) & (qc_ht.f_stat <= 0.6) & (hl.is_defined(
                               qc_ht.normalized_y_coverage) & (qc_ht.normalized_y_coverage > 0.1))),
                           sex_aneuploidy=(qc_ht.f_stat < 0.4) & hl.is_defined(qc_ht.normalized_y_coverage) & (qc_ht.normalized_y_coverage > 0.1))

    sex_expr = (hl.case()
                .when(qc_ht.ambiguous_sex, "ambiguous_sex")
                .when(qc_ht.sex_aneuploidy, "sex_aneuploidy")
                .when(qc_ht.is_female, "female")
                .default("male"))

    qc_ht = qc_ht.annotate(
        sex=sex_expr, data_type='exomes').key_by('data_type', 's')

    qc_ht_file = mtdir + "mt_ambiguous_sex_samples.ht"
    qc_ht.write(qc_ht_file, overwrite=True)


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

    mt_in_file = mtdir + "/gatk_unprocessed.mt"
    print("Reading input matrix")
    mt_unfiiltered = hl.read_matrix_table(mt_in_file)

    #apply hard fitlers
    mt_filtered = apply_hard_filters(mt_unfiiltered, mtdir)

    #impute sex and annotate samples wih ambiguous sex
    mt_sex = impute_sex(mt_filtered, mtdir, annotdir, male_threshold=0.6)
    # annotate_ambiguous_sex(mt_sex, mtdir)


if __name__ == '__main__':
    main() 
