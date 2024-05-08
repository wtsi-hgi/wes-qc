#!/usr/bin/env python
"""
This is the wrapper Python script executing Snakemake pipeline for the variant QC
We need this wrapper script to run Snakemake inside Spark.

Usage: hlrun_local wesqc.py args
"""
import re
import sys
from snakemake import main

if __name__ == "__main__":
    sys.argv[0] = re.sub(r"(-script\.pyw|\.exe)?$", "", sys.argv[0])
    add_args = ["-p", "--cores", "1"]
    if not "-s" in sys.argv:
        add_args += ["-s", "wesqc.Snakefile"]
    sys.argv = [sys.argv[0]] + add_args + sys.argv[1:]
    sys.exit(main())
