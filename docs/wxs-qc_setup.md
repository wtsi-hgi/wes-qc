# How to set up and run the WxS-QC pipeline

The pipeline code is infrastructure agnostic and can run on a single machine.
However, Hail data structures are “heavy” in terms of storage space and memory.
Before starting, ensure you have enough memory and storage space to proceed.

If you're working inside the Wellcome Sanger Institute infrastructure,
we recommend that you use `osdataproc` to create cluster with Hail and Spark backend.
Follow the [Hail on SPARK](https://hgi-projects.pages.internal.sanger.ac.uk/documentation/docs/tutorials/hail-on-spark/)
guide to create such a cluster.


## Set up the codebase

Clone the repository using:
```shell
git clone https://github.com/wtsi-hgi/wxs-qc.git
cd wxs-qc
```

### Set up the environment (local installation only)

If you are running the code on a local machine,
set up virtual environment using `uv`.

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh # Install `uv` system-wide
```

Create and activate your virtual environment
```bash
uv venv
source .venv/bin/activate
uv sync
```

Set up environment variables:
add the project folder to the `PYTHONPATH`,
set memory size for SPARK.

```bash
export PYTHONPATH=$( pwd ):$PYTHONPATH
export PYSPARK_SUBMIT_ARGS="--driver-memory 48G --executor-memory 48G pyspark-shell"
```

**Note:** You need to activate venv and modify variables each time you open a new shell.

## How to run the code

### Manually running the code on a local machine

To manually run the code on a local machine,
activate the environment, run the Python, and provide the path to the pipeline script:

```shell
python 1-import_data/1-import_gatk_vcfs_to_hail.py
```

### Manually running the code on a Hail cluster

To submit the jobs on a Hail cluster, you need to set up environment variables
to include the directory you originally cloned the git repo into.
To get a correct Python path, you need to have the virtual environment activated.

```shell
export PYTHONPATH=$PYTHONPATH:$(pwd)
export PYSPARK_PYTHON=$(which python)
export PYSPARK_DRIVER_PYTHON=$(which python)
```

No you can run the pipeline script via `spark-submit`:

```shell
spark-submit 1-import_data/1-import_gatk_vcfs_to_hail.py
```

We suggest running all code on a cluster in `tmux`/`screen` session to avoid
script termination in case of any network issues.


### Automatically syncing and running the code on the cluster from a local machine

If you want to modify the code on your local machine,
and then run it on the cluster, you can use two scripts provided
in the `scripts` folder.

* `hlrun_local` - runs the Python script via `spark-submit`.
* `hlrun_remote` - runs the code on the Spark cluster from your local machine.
  It performs a series of operations:
  * Syncs the codebase to the remote cluster, defined by the environment variable `$hail_cluster`.
    The variable can contain the full host definition (`user@hostname`) or only hostname from the SSH config file.
  * Creates `tmux` session on the remote cluster.
  * Runs the Python script via `hlrun_local`
  * Attaches to the tmux session to monitor the progress

**Warning**

The `hlrun_remote` is designed to work with only one tmux session.
To start a new task via `hlrun_remote`, first end the existing tmux session, if it exists.


### Running the code via Jupyter notebook

You can run the code in the provided Jupyter notebook
where all the steps are arranged in a sequence and divided into sections
(e.g. 0-resource_preparation, 1-import_data, 2-sample_qc, 3-variant_qc, 4-genotype_qc).

The notebook is located as `scripts/run-wxs-qc-pipeline-all-steps.ipynb`.

It uses `hlrun_local` to run the code, which will output the log file to the current directory,
with the prefix of the step name, e.g. `hlrun_3-1-generate-truth-sets_20250102_125729.log`.

For details, refer to the Markdown comments in the notebook.

