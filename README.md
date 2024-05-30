# Sample and variant QC for exome cohort study

Based on the gnomad hail-based QC pipelines and on Pavlos Antoniou's code
in https://github.com/wtsi-hgi/exomeQC-hail-gnomad/tree/lustre

The current codebase has several branches - one for each dataset.
The unification of branches and the code refactoring is in progress.

The howto for the code is here:
https://hgi-projects.pages.internal.sanger.ac.uk/documentation/docs/how-to-guides/wes-qc-hail/

# How to run the code

The code is designed to run on the SPARK cluster with the
[Hail library](https://hail.is/) installed.
The manula for the cluster setting up is
[here](https://hgi-projects.pages.internal.sanger.ac.uk/documentation/docs/tutorials/hail-on-spark/#destructing-and-re-creating-cluster).

To run code, two scripts are used:
* `hlrun_local` - runs the Python script via `spark-submit` on the cluster.
* `hlrun_remote` - a wrapper to run the code on the cluster directly from local developer's machine
  It performs a set of operation:
  * Sync the codebase to the remote cluster
  * Create **tmux** session on the remoter cluster
  * Run the Python script in the **tmux** session via `hlrun_local`
  * Attach to the **tmux** session to monitor the progress

**Warning**
The `hlrun_remote` script is designed to work with only one tmux session.
To start a new task via `hlrun_remote`, first exit the existing tmux session, if it exists.

Example:

```shell
# Ensure that you're in the main repo folder
hlrun_remote 4-genotype_qc/1a-apply_range_of_hard_filters.py
```

The script takes the address of Hail cluster form the environment variable `$hail_cluster`.
You need to set it in your `.bashrc` as follows:

```shell
export hail_cluster="ubuntu@172.27.11.11"
```


## How to run hard filter evaluation

The hard filter evaluation step is a POC of using
[Snakemake](https://snakemake.readthedocs.io/en/stable/#)
pipeline engine as a dependency tracking system.

The main part of pipeline is available in the main folder in the `wesqc.Snakefile`
To run the pipeline you need to start wrapper `wesqc.py`:

```shell
hlrun_remote wesqc.py
```

To see the set of tasks to execute, run the pipeline with `-n` switch.

To understand the basics of Snakemake workflow management (it's quite different from "regular" pipelines),
see the [Snakemake short intro](https://drive.google.com/file/d/1QgJmc0V0zIU5ydmft3RMyPFQsgqXyY6N/view?usp=sharing).

For all questions, please contact [Gennadii Zakharov](mailto:gz3@sanger.ac.uk)
