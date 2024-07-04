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


def change_ids_in_csv(input_csv: str, map_csv: str, output_csv: str) -> None:
    # Load the mapping file into a dictionary
    id_map = {}
    with open(map_csv, mode="r") as map_file:
        map_reader = csv.DictReader(map_file)
        for row in map_reader:
            id_map[row["OldID"]] = row["NewID"]

    # Open the input CSV file and the output CSV file
    with open(input_csv, mode="r") as infile, open(output_csv, mode="w", newline="") as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames if reader.fieldnames is not None else []
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)

        # Write the header to the output CSV file
        writer.writeheader()

        # Process each row in the input CSV file
        for row in reader:
            old_id = row["SampleID"]
            if old_id in id_map:
                row["SampleID"] = id_map[old_id]
            else:
                print(f"!!! WARNING: No new ID found for {old_id} in the map file.")
            writer.writerow(row)


def main() -> None:
    inputs = parse_config()

    data_root: str = inputs["data_root"]

    annot_dir: str = os.path.join(data_root, inputs["annotation_lustre_dir"])

    fam_files: list[str] = [os.path.join(annot_dir, fam_file) for fam_file in inputs["fam_pedigree"]]
    print("=== Extracting sex information from pedigree files: ===")
    print("\n".join(fam_files))

    self_reported_sex: str = os.path.join(annot_dir, inputs["self_reported_sex"])
    extract_sex_from_multiple_fam(fam_files, self_reported_sex)

    if inputs["fam_id_remap"] is not None:
        temp_self_reported_sex = os.path.join(annot_dir, "temp_self_reported_sex.csv")
        id_map = os.path.join(annot_dir, inputs["fam_id_remap"])

        os.rename(self_reported_sex, temp_self_reported_sex)
        change_ids_in_csv(temp_self_reported_sex, id_map, self_reported_sex)
    print(f"Saving sex annotation to {self_reported_sex}")


if __name__ == "__main__":
    main()
