#apply gnomad's hard filters and impute sex
#input gatk_unprocessed.mt from step 1.1
import hail as hl
import pyspark
from utils.utils import parse_config


def apply_hard_filters(mt: hl.MatrixTable, mtdir: str) -> str:
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

    return filtered_mt_file


def impute_sex(filename: str, mtdir: str, annotdir: str, male_threshold: float = 0.8, female_threshold: float = 0.5) -> str:
    '''
    Imputes sex, exports data, and annotates mt with this data
    :param MatrixTable mt: MT containing samples to be ascertained for sex
    :param str mtdir: directory output matrix tables are written to
    :param str annotdir: directory annotation files are written to
    :return: MatrixTable with imputed sex annotations stashed in column annotation 'sex_check'
    :rtype: MatrixTable
    '''
    print("Imputing sex with male_threshold = " + str(male_threshold) + " and female threshold = " + str(female_threshold))
    mt = hl.read_matrix_table(filename)

    #filter to X and select unphased diploid genotypes - no need to filter to X as impute_sex will do this
    #mt1 = hl.filter_intervals(mt, [hl.parse_locus_interval('chrX')])
    mt1 = hl.split_multi_hts(mt)
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

    return sex_mt_file


def identify_inconsistencies(filename: str, mtdir: str, annotdir: str, resourcedir: str, metadata_file: str):
    '''
    Find samples where annotated sex conflicts with the sex in our metadata
    Find samples where sex is not annotated
    Find samples where f_stat is between 0.2 and 0.8
    :param MatrixTable mt: MT containing imputed sex in column 'sex_check'
    :param str mtdir: directory output matrix tables are written to
    :param str annotdir: directory annotation files are written to
    :param str resourcedir: directory annotation files are written to
    '''
    print("Annotating samples with inconsistencies:")
    mt = hl.read_matrix_table(filename)
    qc_ht = mt.cols()
    #convert is_female boolean to sex
    sex_expr = (hl.case()
            .when(qc_ht.is_female, "female")
            .when(qc_ht.is_female == False, "male")
            .default("undetermined"))
    qc_ht = qc_ht.annotate(sex=sex_expr).key_by('s')

    #annotate with manifest sex - keyed on ega to match identifiers in matrixtable
    metadata_ht = hl.import_table(metadata_file, delimiter="\t", no_header=False)

    #we only want those from the metadata file where sex is known
    metadata_ht = metadata_ht.filter((metadata_ht.sex == '1') | (metadata_ht.sex == '2'))

    #annotate the sex-predictions with the manifest sex annotation - need to use a join here
    metadata_ht = metadata_ht.key_by('sample_id')
    ht_joined = qc_ht.annotate(manifest_sex = metadata_ht[qc_ht.s].sex)

    #identify samples where imputed sex and manifest sex conflict
    conflicting_sex_ht = ht_joined.filter(((ht_joined.sex == 'male') & (ht_joined.manifest_sex == '2')) | (
        (ht_joined.sex == 'female') & (ht_joined.manifest_sex == '1')))
    conflicting_sex_ht.export(annotdir + '/conflicting_sex.txt')

    #identify samples where f stat is between 0.2 and 0.8
    f_stat_ht = qc_ht.filter( (qc_ht.f_stat > 0.2) & (qc_ht.f_stat < 0.8) )
    f_stat_ht.export(annotdir + '/sex_annotation_f_stat_outliers.txt.bgz')
    

def main():
    #set up
    inputs = parse_config()
    importmtdir = inputs['load_matrixtables_lustre_dir']
    mtdir = inputs['matrixtables_lustre_dir']
    annotdir = inputs['annotation_lustre_dir']
    resourcedir = inputs['resource_dir']

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

    #impute sex
    mt_sex = impute_sex(mt_filtered, mtdir, annotdir, male_threshold=0.79, female_threshold=0.55)

    # annotate_ambiguous_sex(mt_sex, mtdir)
    metadata_file = inputs['metadata_file']
    #the file has two columns 'sample_id' and 'sex'
    identify_inconsistencies(mt_sex, mtdir, annotdir, resourcedir, metadata_file)


if __name__ == '__main__':
    main() 
