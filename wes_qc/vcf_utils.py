import os.path
import re
from typing import List, Dict

chr_name_pattern = r"^chr([0-9XY][0-9]?)_([0-9]+)-([0-9]+)"

chr_to_alphabet: Dict[str, str] = {"X": "5X", "Y": "6Y", "M": "9M"}


def sort_vcf_files(vcf_files: List[str]) -> List[str]:
    """ """

    def extract_key(filename: str) -> tuple[str, int]:
        # Extract chromosome and position
        match = re.match(chr_name_pattern, os.path.basename(filename))
        if match:
            chrom = match.group(1)
            pos = int(match.group(2))
            # Convert chromosomes to str leading zero to allow correct order
            if chrom not in ("X", "Y", "M"):
                chrom = f"{int(chrom):02d}"
            else:
                chrom = "3" + chrom
            return (chrom, pos)
        else:
            raise ValueError(f"Can't sort file name: {filename}")

    # Sort files by extracted key
    sorted_files = sorted(vcf_files, key=extract_key)
    return sorted_files


if __name__ == "__main__":
    vcf_files = [
        "chrY_19744285-19744634.sort.pre-vep.vcf.gz",
        "chrY_284091-284410.sort.pre-vep.vcf.gz",
        # "chrM_1-16569.sort.noaltcontig.vcf.gz",
        "chrX_91876681-91879373.sort.pre-vep.vcf.gz",
        "chrX_9745939-9746258.sort.pre-vep.vcf.gz",
        "chr1_12345-67890.sort.pre-vep.vcf.gz",
        # "chr1_KI270711v1_random_9265-9585.sort.vcf.gz",  # Chromosome with altcontig
        "chr2_12345-67890.sort.pre-vep.vcf.gz",
        "chr20_2345-67890.sort.pre-vep.vcf.gz",
        "chr10_12345-67890.sort.pre-vep.vcf.gz",
    ]
    print("\n".join(sort_vcf_files(vcf_files)))
