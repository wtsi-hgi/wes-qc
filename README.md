# Sample and variant QC for exome cohort studies

This repository contains the code
used for the QC of human exome cohorts

The code is written by members of **Wellcome Sanger HGI group**
(https://www.sanger.ac.uk/group/human-genetics-informatics-hgi/)
based on the **gnomAD QC** pipeline by **Broad institute**
(https://github.com/broadinstitute/gnomad_qc/tree/main).

The current codebase has several branches — one for each dataset.
The unification of branches and the code refactoring are in progress.

The howto for the code is here:
https://hgi-projects.pages.internal.sanger.ac.uk/documentation/docs/how-to-guides/wes-qc-hail/

# How to run the code and QC your data

The brief howto for the QC process is available in the
[WES-QC Hail howto](docs/wes-qc-hail.md).

The code is designed to run on the SPARK cluster with the
[Hail library](https://hail.is/) installed.
The manual for the cluster setting up is
[here](docs/hail-on-spark.md).

# How to get the latest changes from main

When working on an analysis branch, you can retrieve the latest changes from main by running:
```bash
make update
```
This will fetch the latest changes from main and rebase the current branch onto it.

# Developer's howto

## How to run the tests and calculate coverage

The tests currently require running on the SPARK cluster. There are plans to make them runnable locally.
They can be run by commands defined in `Makefile`.

To run all the tests:
```bash
make test
```
Or you can specify the type of test to run
```bash
make unit-test
make integration-test
```

To run the tests with coverage:
```bash
make unit-test-coverage
make integration-test-coverage
```

## To run pre-commit hooks on commit

1. Install pre-commit
```shell
pip install pre-commit
```
2. `pre-commit` will automatically run on every commit
3. To run pre-commit manually on specific files
```shell
pre-commit run --files <file1> <file2>
```
4. `mypy` is configured to run manually because now it produces too many errors. To run it:
```shell
pre-commit run --hook-stage manual
```
