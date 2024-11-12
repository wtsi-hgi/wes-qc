# Sample and variant QC for exome cohort studies

This repository contains the code
used for the QC of human exome cohorts

The code is written by members of **Wellcome Sanger HGI group**
(https://www.sanger.ac.uk/group/human-genetics-informatics-hgi/)
based on the **gnomAD QC** pipeline by **Broad institute**
(https://github.com/broadinstitute/gnomad_qc/tree/main).

The current codebase has several branches - one  for each dataset.
The unification of branches and the code refactoring is in progress.

The howto for the code is here:
https://hgi-projects.pages.internal.sanger.ac.uk/documentation/docs/how-to-guides/wes-qc-hail/

# How to run the code

The code is designed to run on the SPARK cluster with the
[Hail library](https://hail.is/) installed.
The manual for the cluster setting up is
[here](https://hgi-projects.pages.internal.sanger.ac.uk/documentation/docs/tutorials/hail-on-spark/#destructing-and-re-creating-cluster).

To run the code on a Spark cluster, two scripts are used:
* `hlrun_local` - runs the Python script via `spark-submit`. You need to run it on the spark master node on your cluster.
* `hlrun_remote` - runs the code on the Spark cluster form your local machine.
It performs a series of operations:
  * Sync the codebase to the remote cluster
  * Create tmux session on the remoter cluster
  * Run the Python script via `hlrun_local`
  * Attach to the tmux session to monitor the progress

**Warning**

The `hlrun_remote` is designed to work with only one tmux session.
To start a new task via `hlrun_remote`, first end the existing tmux session, if it exists.


# How to run the tests and with coverage

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

# Developer's howto

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
4. `mypy` is configured to run manually because now it produce too many errors. To run it:
```shell
pre-commit run --hook-stage manual
```
