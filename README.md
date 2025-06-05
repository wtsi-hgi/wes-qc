# WxS-QC - sample and variant QC for genome/exome cohort studies

This repository contains the pipeline
for quality control of human genome/exome cohorts

The code is written by members of **Wellcome Sanger HGI group**
(https://www.sanger.ac.uk/group/human-genetics-informatics-hgi/)
based on the **gnomAD QC v3** pipeline by **Broad institute**
(https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v3).

The code is designed to run standalone or on the SPARK cluster with the
[Hail library](https://hail.is/) installed.

## How to QC your data

The high-level description of the QC process is available in a separate document: 
[WxS-QC concepts](docs/wxs-qc_concepts.md).

Detailed howto for the QC process is here:
[WcS-QC howto](docs/wxs-qc_howto.md).

### How to get the latest changes from the `main` branch

When working on an analysis branch, you can retrieve the latest changes from the `main` branch by running:

```bash
make update
```

This will fetch the latest changes from the `main` branch and rebase your current branch onto it.
If there are any unstaged changes in the branch, you will be asked to commit or stash them first.

## Developer's howto

For improving the pipeline and developing new functionality,
the [Developer's howto](docs/wxs-qc_development.md) is available.