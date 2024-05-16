# Sample and variant QC for exome cohort study

Based on the gnomad hail-based QC pipelines and on Pavlos Antoniou's code
in https://github.com/wtsi-hgi/exomeQC-hail-gnomad/tree/lustre

The current codebase has several branches - one  for each dataset.
The unification of branches and the code refactoring is in progress.

The howto for the code is here:
https://hgi-projects.pages.internal.sanger.ac.uk/documentation/docs/how-to-guides/wes-qc-hail/

# How to run the code

The code is designed to run on the SPARK cluster with the
[Hail library](https://hail.is/) installed.
The manula for the cluster setting up is
[here](https://hgi-projects.pages.internal.sanger.ac.uk/documentation/docs/tutorials/hail-on-spark/#destructing-and-re-creating-cluster).

To run code, two scripts are used:
* `hlrun-local` - runs the Python script via `spark-submit`.
* `hlrun-remote` - makes a series of operation:
  * Sync the codebase to the remote cluster
  * Create tmux session on the remoter cluster
  * Run the Python script via `hlrun-local`
  * Attach to the tmux session to monitor the progress

**Warning**
The hlrun-remote designed to work with only one tmux session.
To start a new task via `hlrun-remote`, first exit the existing tmux session, if it exists.
