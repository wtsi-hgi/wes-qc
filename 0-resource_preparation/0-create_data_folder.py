import os

from utils.utils import parse_config


def main() -> None:
    # = STEP SETUP = #
    config = parse_config()
    conf = config["general"]

    # = STEP PARAMETERS = #
    data_root = conf["data_root"]
    metadata_dir = conf["metadata_dir"]
    annotation_dir = conf["annotation_dir"]
    matrixtables_dir = conf["matrixtables_dir"]
    resource_dir = conf["resource_dir"]
    onekg_resource_dir = conf["onekg_resource_dir"]
    plots_dir = conf["plots_dir"]
    var_qc_rf_dir = conf["var_qc_rf_dir"]
    vcf_export_dir = conf["vcf_export_dir"]
    gatk_vcf_indir = config["step1"]["gatk_vcf_indir"]

    # = STEP LOGIC = #
    os.makedirs(data_root, exist_ok=True)
    for directory in [
        metadata_dir,
        annotation_dir,
        matrixtables_dir,
        resource_dir,
        onekg_resource_dir,
        plots_dir,
        var_qc_rf_dir,
        vcf_export_dir,
        gatk_vcf_indir,
    ]:
        os.makedirs(directory, exist_ok=True)


if __name__ == "__main__":
    main()
