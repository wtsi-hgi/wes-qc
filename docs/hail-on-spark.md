# Getting Started With Hail On Openstack

This tutorial describes how to initiate a cluster on **Openstack** suitable for running **Hail**.

The command line tool [osdataproc](https://github.com/wtsi-hgi/osdataproc/) is used to create an openstack cluster with [Apache Spark](https://spark.apache.org/)  and [Apache Hadoop](https://hadoop.apache.org/) configured. The cluster is provisioned with [Jupyterlab](https://jupyter.org/) for interactive computing, [Hail](https://hail.is/) for genomic analysis, and [Netdata](https://www.netdata.cloud/) for monitoring.
## Install `osdataproc`
This part of work is done on your local machine.

> [!NOTE]
The `osdataproc` utility requires python 3.9. Before proceeding, ensure that you're using the correct Python version. Ensure that the Python version is correct:
>```shell
>python --version
>```
>
>The easiest way to install the correct Python version is using [Conda](https://docs.conda.io/en/latest/).
>```shell
>conda create -n py39 python=3.9
>conda activate py39
> ```
>or  install the correct Python version locally:
>``` shell
>wget https://www.python.org/ftp/python/3.9.8/Python-3.9.8.tgz
>tar -zxvf Python-3.9.8.tgz
>cd Python-3.9.8/
>mkdir ~/.localpython
># Prepare the environment for building
>./configure --prefix=$HOME/.localpython
># Building the system
> make
> # Implement the installation
>make install
>```
>Running the code above will install python 3.9.8 at `~/localpython/bin/python3.9`

### Create and activate the Python virtual environment for `osdataproc`

Set up a virtual environment with newly installed Python.
```bash
virtualenv osdataproc_env -p path/to/python3.9
# activate
source osdataproc_env/bin/activate
# deactivate with simply `deactivate`
```

### Install Terraform

`osdataproc` requires [Terraform](https://www.terraform.io/) version 1.4 or higher. Find the location of the latest Terraform [here](https://developer.hashicorp.com/terraform/downloads?product_intent=terraform), download and unzip into a location on your path. The easy way is to use the `bin` folder of the environment created above. You should download archive for amd64, it works on M1/M3 MacBooks.

```bash
cd osdataproc_env/bin
wget https://releases.hashicorp.com/terraform/1.9.3/terraform_1.9.3_linux_386.zip
unzip terraform_1.8.2_darwin_amd64.zip
cd ../..
```

### Install osdataproc

Clone the `osdataproc` git repository in a separate folder on your machine, and install it via `pip`.

```bash
git clone https://github.com/wtsi-hgi/osdataproc.git
cd osdataproc
pip install -e .
cd ..
```

Check, that `osdataproc` installed successfully:
```shell
osdataproc create --help
```

> [!WARNING]
Terraform stores all cluster configuration data in the  `./terraform/terraform.tfstate.d` folder in the `osdataproc` folder. **Don't remove it**. Otherwise, you lose access to cluster configuration and won't be able to suspen/resume/destroy created clusters.

>[!TIP]
>If when installing `osdataproc` the wheel fail to be built, try updating `wheel` and `setuptools` packages.
>``` shell
>pip install wheel -U # 0.43.0 version
>pip install setuptools -U # 71.1.0 version
>```

## Managing a cluster
### Create a cluster

For naming your clusters and volumes, we recommend prefixing with your username, similarly for your volumes (nfs-volume).
**Example:** `gz3-hail`

You will be prompted for a password, this password is used for accessing Jupyter and Spark master node. Since you may have to swap clusters or allow other people to access this password needs to be something you are happy to share with others.

**Example**: To create a cluster using 50 `m2.medium` workers
and create a new volume called `gz3-hail` run the following script:

```shell
eval `ssh-agent -s` \
ssh-add path/to/public/key

osdataproc create [--num-workers]    <Number of desired worker nodes>
                  [--public-key]     <Path to public key file>
                  [--flavour]        <OpenStack flavour to use>
                  [--network-name]   <OpenStack network to use>
                  [--lustre-network] <OpenStack Lustre provider network to use>
                  [--image-name]     <OpenStack image to use - Ubuntu images only>
                  [--nfs-volume]     <Name/ID of volume to attach or create as NFS shared volume>
                  [--volume-size]    <Size of OpenStack volume to create>
                  [--device-name]    <Device mountpoint name of volume>
                  [--floating-ip]    <OpenStack floating IP to associate to master node - will automatically create one if not specified>
                  <cluster_name>
```

Wait until Terraform completes cluster creation. The creation of clusters can take some time and needs your laptop to be connected to the network.

>[!IMPORTANT]
>The IP address of your new cluster will be printed to STDOUT. **Save it somewhere**.

When the cluster completion finishes, inspect the log. There should be no failed tasks:
### Destructing and re-creating cluster
There are many reasons, why the cluster creation can fail.
In most cases, it's easier to delete and re-create the cluster.

To delete cluster, run the following command:
```bash
osdataproc destroy "${cluster_name}"
```

To re-create cluster, run the same command without specifying `--volume-size`. The new cluster will use existing volume. In most cases, volumes are created correctly, so you don't need to delete and re-create it.
## Running Hail tasks on the cluster

Log into your new cluster using the IP address created when your new cluster was created. For example:

```bash
ssh ubuntu@<public_ip>
```

### Installing additional software

The WES-QC code requires Python GnomaAD package. To install it, run the following:

```shell
sudo apt install postgresql python3.9-dev libpq-dev
pip install gnomad
```

### Check cluster status
To check the cluster status use the SPARK web interface:
<http://<public_ip>/spark/>

### Submit Hail scripts

Hail scripts are submitted to worker nodes using spark as follows:

```bash
export PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python
spark-submit /path/to/hail_script.py
```

Alternatively, you can activate environment and use spark commands directly

```bash
source /home/ubuntu/venv/bin/activate
spark-submit /path/to/hail_script.py
```

It is highly advisable to open `tmux` or `screen` session before submitting spark jobs.

### Access Jupyter

Jupyter can be accessed via a web browser as follows (using the IP address created when you created your new cluster):
<https://<public_ip>/jupyter/lab>

You will be prompted for a password that you used to create the cluster.

>[!WARNING]
All Hail tasks (both command-line and interactive via Jupyter) occupy all working nodes. You can't run Jupyter notebook and command-lite script simultaneously.
>
>To kill a Jupyter job, shut down all Jupyter kernels, or kill the `Hail` process via Spark master web interface.

<!--### Workaround netdata log issue
Often the netdata access logs grow too fast. To avoid running out of free space, we suggest changing logrotate settings.

First, change the logrotate schedule:
```shell
sudo ln -s /etc/cron.daily/logrotate /etc/cron.hourly/logrotate
```

Then, add the following lines to the `/etc/logrotate.d/netdata` file:
```
/var/log/netdata/*.log {
        maxsize 1G
        hourly
        rotate 1
        ...
}
```

Finally, run Logrotate manually to ensure that it works without errors
```shell
sudo logrotate /etc/logrotate.d/netdata
```

## Troubleshooting

### Cleaning up Netdata logs

In case if you cluster reports `No space left on the device`
the most probable source of this issue is netdata logs.
To inspect it, go to the netdata folder, and check the log files size:
```shell
cd /var/log/netdata
ls -lh
```

If the size of netdata logs is too big, remove it:
```shell
sudo rm -rf *.log.*
```

### Resetting HDFS safe mode

In case of an issue with the free space on the device (both caused by local FS overflow or by Lustre quota)
the HDFS will switch in the safe mode. The Hail processing will throw the error message:
`Caused by: org.apache.hadoop.ipc.RemoteException(org.apache.hadoop.hdfs.server.namenode.SafeModeException):
Cannot create directory /shared/spark-logs. Name node is in safe mode.
`

To continue working with Hail, you need to manually move it bach to operational mode. To do it:
1. Deal with the absence of free space (clean up logs, temporary folders, unused matrixtables, etc.).
2. Turn HDFS bach to the operational mode:

```shell
hdfs dfsadmin -safemode leave
```

### Manual cluster cleanup
In the case of cluster creation/destruction failure (for example, due to connection loss),
`osdataproc` may not be able to clean up all cluster resources.

To manually clean up all requested resources via Theta, do the following, do the following:

>[!WARNING]
Working under the project user, you can delete resources belonging to other project members.
>
>Please, be cautious and double-check every step.

1. Log in to [Theta](https://theta.internal.sanger.ac.uk/project/)
2. Choose the project under `Projects`
3. Clean up instances
   1. Navigate to `instances`
   2. Set the filter by the instance name and find all instances from your cluster
   3. **Double check, that you selected only instances belonging to your cluster**
   4. Remove instances
4. Clean up Lustre network ports
   1. Navigate to `Networks`
   2. Choose the secure Lustre network
   3. Go to the `Ports` tab
   4. Set the filter and find all ports from your cluster
   5. **Double check, that you selected only ports belonging to your cluster**
   6. Remove ports
5. Clean up Cloudforms network
   1. Navigate to `Networks`
   2. Choose the Cloudforms network
   3. Remove ports the same way as yu did in point 4
6. Remove cluster security groups:
   1. Navigate to `Security groups`
   2. Set the filter two security groups created for the cluster.
   3. Remove security groups -->
