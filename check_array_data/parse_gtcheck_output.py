#!/usr/bin/env python3
# parse output of bcftools gtcheck
import gzip
import argparse
from typing import Any

MatchDict = dict[str, dict[str, Any]]


def parse_gtcheck_file(gtcheck_file: str) -> MatchDict:
    """
    parse gtcheck file and identify best hit
    """
    best_hits: MatchDict = {}
    if gtcheck_file.endswith(".gz"):
        file_context = gzip.open(gtcheck_file, "rt")
    else:
        file_context = open(gtcheck_file)
    with file_context as f:
        for line in f:
            if not line.startswith("DC"):
                continue
            linedata = line.split()
            sample = linedata[1]
            famid = linedata[2]
            discordance = float(linedata[3])
            nsites = int(linedata[5])
            if nsites > 0:
                score = discordance / nsites
            else:
                score = 1

            if sample in best_hits.keys():
                if score < best_hits[sample]["score"]:
                    best_hits[sample] = {"fam": famid, "discordance": discordance, "nsites": nsites, "score": score}
            else:
                best_hits[sample] = {"fam": famid, "discordance": discordance, "nsites": nsites, "score": score}

    return best_hits


def write_output(outfile: str, best_hits: MatchDict) -> None:
    """
    write output file
    """
    print("Writing output")
    with open(outfile, "w") as o:
        header = ("\t").join(["#EGA", "broad_id_top_hit", "sum_discordance", "nsites", "discordance/nsites"])
        o.write(header + "\n")
        for sample in best_hits.keys():
            outline = ("\t").join(
                [
                    sample,
                    best_hits[sample]["fam"],
                    str(best_hits[sample]["discordance"]),
                    str(best_hits[sample]["nsites"]),
                    str(best_hits[sample]["score"]),
                ]
            )
            o.write(outline + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("gtcheck_file", help="gtcheck output file")
    parser.add_argument("output_file", nargs="?", help="output file")
    args = parser.parse_args()
    if not args.output_file:
        args.output_file = args.gtcheck_file + ".parsed.txt"

    best_hits = parse_gtcheck_file(args.gtcheck_file)
    write_output(args.output_file, best_hits)


if __name__ == "__main__":
    main()
