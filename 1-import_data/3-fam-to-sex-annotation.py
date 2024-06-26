import csv
import os
from typing import Dict

from utils.utils import parse_config


def extract_sex_from_multiple_fam(fam_files: list[str], output_csv: str) -> None:
    # Define the sex mapping based on PLINK codes
    sex_map: Dict[str, str] = {"1": "Male", "2": "Female", "0": "Unknown"}

    # Set to track duplicate sample IDs
    seen_sample_ids = set()

    with open(output_csv, "w", newline="") as out_csv:
        csv_writer = csv.writer(out_csv)

        # Write the header for the output CSV
        csv_writer.writerow(["SampleID", "Sex"])

        for fam_file in fam_files:
            with open(fam_file, "r") as fam:
                fam_reader = csv.reader(fam, delimiter=" ")

                for row in fam_reader:
                    sample_id = row[1]
                    sex_code = row[4]
                    sex = sex_map.get(sex_code, "Unknown")

                    if sample_id in seen_sample_ids:
                        raise ValueError(f"Warning: Duplicate sample ID found: {sample_id}")
                    else:
                        seen_sample_ids.add(sample_id)
                        csv_writer.writerow([sample_id, sex])


def main() -> None:
    inputs = parse_config()

    data_root: str = inputs["data_root"]

    annot_dir: str = os.path.join(data_root, inputs["annotation_lustre_dir"])

    fam_files: list[str] = [os.path.join(annot_dir, fam_file) for fam_file in inputs["fam_pedigree"]]
    output_csv: str = os.path.join(annot_dir, inputs["sex_and_ethnicity_annotation_file"])
    print("\n".join(fam_files))

    extract_sex_from_multiple_fam(fam_files, output_csv)


if __name__ == "__main__":
    main()
