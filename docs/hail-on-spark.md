---
tags:
  - Tutorial
  - Compute
  - OpenStack
  - Hail
---

# Getting Started With Hail On Openstack

!!! warning
    This documentation is under development and may be incomplete.

This tutorial describes how to initiate a cluster on **Openstack** suitable for running **Hail**. 

The command line tool [osdataproc](https://github.com/wtsi-hgi/osdataproc/) 
is used to create an openstack cluster with [Apache Spark](https://spark.apache.org/) 
and [Apache Hadoop](https://hadoop.apache.org/) configured. 
The cluster is provisioned with [Jupyterlab](https://jupyter.org/) for interactive computing, 
[Hail](https://hail.is/) for genomic analysis, and [Netdata](https://www.netdata.cloud/) for monitoring. 

## Before you start

Before you start, ensure that you have created an ssh keypair with [ssh-keygen](https://manpages.ubuntu.com/manpages/xenial/man1/ssh-keygen.1.html).
Use the RSA key with blank password. See the 
[SSG confluence page](https://ssg-confluence.internal.sanger.ac.uk/display/FARM/All+things+SSH)
for details.

```shell
ssh-keygen -t rsa
```

You need access to an openstack project in order to create a cluster.
To obtain your OpenStack password, follow the 
[manual](https://ssg-confluence.internal.sanger.ac.uk/pages/viewpage.action?pageId=66031299).
You don't need to reset the password, you can copy and use the auto-generated one.

Upload your RSA ssh key to OpenStack following the 
[manual](https://ssg-confluence.internal.sanger.ac.uk/display/OPENSTACK/Creating+your+first+machine%2C+using+manageIQ).
You don't need to create a virtual machine, just upload your SSH key.

## Install `osdataproc`
This part of work is done on your local machine.

!!! warning
    The `osdataproc` utility requires python 3.9. 
    Before proceeding, ensure that you're using the correct Python version.
    The easiest way to install the correct Python version is using [Conda](https://docs.conda.io/en/latest/).
    ```shell
    conda create -n py39 python=3.9
    conda activate py39
    ```

Ensure that the Python version is correct:
```shell
python --version                                                                                                                                                                [14:43:21]
Python 3.9.19
```

### Create and activate the Python virtual environment for `osdataproc`

```bash
python3 -m venv osdataproc_env
source osdataproc_env/bin/activate
```

### Install Terraform

`osdataproc` requires [Terraform](https://www.terraform.io/) version 1.4 or higher. 
Find the location of the latest Terraform [here](https://developer.hashicorp.com/terraform/downloads?product_intent=terraform), 
download and unzip into a location on your path.
The easy way is to use the `bin` folder of the environment created above.
You should download archive for amd64, it works on M1/M3 MacBooks.

```bash
cd osdataproc_env/bin
curl https://releases.hashicorp.com/terraform/1.8.2/terraform_1.8.2_darwin_amd64.zip > terraform_1.8.2_darwin_amd64.zip
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

!!! warning
    Terraform stores all cluster configuration data in the `./terraform/terraform.tfstate.d` folder
    in the `osdataproc` folder. **Don't remove it**. 
    Otherwise, you lose access to cluster configuration and won't be able to suspen/resume/destroy 
    created clusters.

## Configure the cluster parameters

### Download and source your project's openrc.sh file

You can find the `openrc.sh` file for your project in [Theta](https://theta.internal.sanger.ac.uk/project/).
Navigate to your project using the drop-down at the top-left part of the page:
![screenshot](/documentation/docs/img/openstack-project.png)
Then download the OpenStack RC file using the drop-down menu at the top right:
![screenshot](/documentation/docs/img/openstack-rc.png)

Source this file
`source /path/to/rcfile.sh`

While you are logged into theta, 
it is a good opportunity to check that there is sufficient resources quote to create your cluster.

### Find the best instance type to use

To choose the instance type, review the available OpenStack instances
[here](https://metrics.internal.sanger.ac.uk/d/000000197/fce-openstack-available-capacity-dynamic?orgId=1&refresh=5m&var-openstack=openstack-theta&var-flavor=All&var-size=All). 
Currently `m1.medium` or `m2.medium` look like good choice.
Hail cluster needs about 40–50 workers.

You can see all available instances on the 
[metrics page](https://metrics.internal.sanger.ac.uk/d/000000197/fce-openstack-available-capacity-dynamic?orgId=1&refresh=5m&var-openstack=openstack-theta&var-flavor=All&var-size=All).
Description of instance types and capabilities is on the 
[SSG confluence](https://ssg-confluence.internal.sanger.ac.uk/display/OPENSTACK/Flavours) page.

### Find a suitable openstack image to use

The only image tested and proved to work with Hail is this one:
`hgi_cant_use_another_image_and_cant_use_an_id_focal-WTSI-lustreonly_200329_c2979847`.
For setting up a Hail cluster, we recommend you to use it. 

Alternatively, you can look for images,   
where the name contains **lustre** 
(has secure lustre capabilities) and **focal** (has python 3.8.10 installed).
To view the available images for your project, use **Compute>Images** menu on Theta:
![screenshot](/documentation/docs/img/openstack-lustre-images.png)

There is documentation about openstack images on 
[SSG confluence](https://ssg-confluence.internal.sanger.ac.uk/display/OPENSTACK/Public+images+-+naming+convention%2C+list%2C+support+and+warnings). 

### Find the name of the secure lustre network

On Theta, navigate to Network>Network Topology,
and view the image of the network to discover the name of the secure lustre network.

## Create a cluster

For naming your clusters and volumes, 
we recommend prefixing with your username, similarly for your volumes (nfs-volume).
Example: `gz3-hail`

You will be prompted for a password, this password is used for accessing Jupyter and Spark master node. 
Since you may have to swap clusters or allow other people to access 
this password needs to be something you are happy to share with others.

**Example**: To create a cluster using 50 m2.medium workers 
and create a new volume called `gz3-hail` run the following script:

```shell
#!/bin/bash

cluster_name="gz3-hail"

osdataproc create           \
  --public-key /Users/gz3/.ssh/id_rsa.pub \
  --num-workers 40                   \
  --flavour m2.medium                \
  --nfs-volume "${cluster_name}-volume"           \
  --volume-size 1000 \
  --image-name hgi_cant_use_another_image_and_cant_use_an_id_focal-WTSI-lustreonly_200329_c2979847 \
  --lustre-network lustre-hgi06      \
 "${cluster_name}" 2>&1 | tee "${cluster_name}.log"
```

Wait until Terraform completes cluster creation.
The creation of clusters can take some time and needs your laptop to be connected to the network.

**The IP address of your new cluster will be printed to STDOUT. Save it somewhere**.

When the cluster completion finishes, inspect the log. There should be no failed tasks:

```shell
PLAY RECAP ***********************************************************************************************
ubuntu@172.27.17.54: ok=83   changed=72   unreachable=0    failed=0    skipped=0    rescued=0    ignored=0 
```

### Destructing and re-creating cluster
There are many reasons, why the cluster creation can fail. 
In most cases, it's easier to delete and re-create the cluster.

To delete cluster, run the following command: `osdataproc destroy "${cluster_name}"`

To re-create cluster, run the same command without specifying `--volume-size`. The new cluster will use existing volume.
In most cases, volumes are created correctly, so you don't need to delete and re-create it.

## Running Hail tasks on the cluster

Log into your new cluster using the IP address created when your new cluster was created. For example:

```bash
ssh ubuntu@172.27.27.98
```

### Installing additional software

The WES-QC code requires Python GnomaAD package. To install it, run the following:

```shell
sudo apt install postgresql python3.9-dev libpq-dev
pip install gnomad
```

### Check cluster status
To check the cluster status use the SPARK web interface:
<http://172.27.17.54/spark/>

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
<https://172.27.27.98/jupyter/lab>

You will be prompted for a password that you used to create the cluster.

!!! warning
    All Hail tasks (both command-line and interactive via Jupyter) occupy all working nodes.
    You can't run Jupyter notebook and command-lite script simultaneously.
    To kill a Jupyter job, shut down all Jupyter kernels, or kill the `Hail` process via Spark master web interface.

### Workaround netdata log issue
Often the netdata access logs grow too fast.
To avoid running out of free space, we suggest changing logrotate settings.

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

!!! warning
    Working under the project user, you can delete resources belonging to other project members.
    
    Please, be cautious and double-check every step.

1. Log in to [Theta](https://theta.internal.sanger.ac.uk/project/)
2. Choose the project under `Projects`
![THETA-project_select](/documentation/docs/img/THETA-project_select.png)
3. Clean up instances
   1. Navigate to `instances`
   2. Set the filter by the instance name and find all instances from your cluster
   3. **Double check, that you selected only instances belonging to your cluster**
   4. Remove instances
![THETA-cleanup-instances](/documentation/docs/img/THETA-cleanup-instances.png)
4. Clean up Lustre network ports
   1. Navigate to `Networks`
   2. Choose the secure Lustre network
   3. Go to the `Ports` tab
   4. Set the filter and find all ports from your cluster
   5. **Double check, that you selected only ports belonging to your cluster**
   6. Remove ports
![THETA-cleanup-ports](/documentation/docs/img/THETA-cleanup-ports.png)
5. Clean up Cloudforms network
   1. Navigate to `Networks`
   2. Choose the Cloudforms network
   3. Remove ports the same way as yu did in point 4
6. Remove cluster security groups:
   1. Navigate to `Security groups`
   2. Set the filter two security groups created for the cluster.
   3. Remove security groups
![THETA-cleanup-security-groups](/documentation/docs/img/THETA-cleanup-security-groups.png)



