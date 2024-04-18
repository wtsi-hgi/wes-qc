#!/usr/bin/env python
import re
import sys
from snakemake import main

if __name__ == '__main__':
    sys.argv[0] = re.sub(r'(-script\.pyw|\.exe)?$', '', sys.argv[0])
    sys.argv = [sys.argv[0], "-p", "--cores", "1"] + sys.argv[1:]
    sys.exit(main())
