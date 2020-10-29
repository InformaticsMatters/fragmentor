# Executing on an AWS ParallelCluster

We can execute on an AWS ParallelCluster environment. A suitable cluster can be
formed by following the instructions in our [nextflow-pcluster] repository.

With a cluster formed you should clone this repository to the **MasterServer**
instance where you can run our playbooks to create a postgres database server
and then execute the fragmentor plays. From here we assume you're
on the cluster's master instance: -

>   You should have your cluster private key in ~/.ssh as described in the
    [nextflow-pcluster] documentation

    $ git clone https://github.com/InformaticsMatters/fragmentor
    $ cd fragmentor/ansible
    $ sudo pip install --upgrade pip
    $ sudo pip install -r ../requirements.txt
    $ ansible-galaxy install -r ../requirements.yaml

Add the cluster private key to the SSH agent: -

    $ eval `ssh-agent`
    $ ssh-add ~/.ssh/nf-pcluster
    
Set your AWS credentials. You'll need these if you're creating the cluster's
postgres database: -

    $ export AWS_DEFAULT_REGION=eu-central-1
    $ export AWS_ACCESS_KEY_ID=?????
    $ export AWS_SECRET_ACCESS_KEY=?????

Life's a lot easier using parameter files with ansible, so create file
that provides variables that satisfy your cluster in order to create
and configure the cluster database. Something like this: -

```yaml
---
db_server_state: present
aws_db_instance_type: t3a.2xlarge
db_volume_size_g: 10
database_cloud_provider: aws
db_shared_buffers_g: 4
db_max_parallel_workers: 8
aws_vpc_subnet_id: <CLUSTER_PUBLIC_SUBNET_ID>
aws_vpc_id: <CLUSTER_VPC_ID>

deployment: production
runpath: /data/share-2/frag
add_backup: no
```

>   Setting `add-backup` to `no` prevents the automatic backup
    from taking place, which typically occurs after the inchi play.

Now create the server: -

    $ ansible-playbook site-db-server.yaml -e @parameters 

Adjust your parameters so that they include the address of the database server.
The server's IP address is printed by the above play: -

    TASK [db-server : Display DB server address (Private IP)] *****************
    Thursday 22 October 2020  18:54:00 +0000 (0:00:00.048)       0:00:24.557 ** 
    ok: [localhost] => {
        "server_result.instances[0].private_ip": "10.0.0.192"
    }

In this case you'd add the following to the parameter file: -

```yaml
database_login_host: 10.0.0.192
```

Using the AWS console wait for the database server instance to become ready
(initialise) before trying to configure it.

You will need to install the `ec2.py` dynamic inventory
script, provided by ansible. The following installs the script locally
and then runs an ansible `ping` to ensure the database server can be found: -

    $ wget https://raw.githubusercontent.com/ansible/ansible/stable-2.9/contrib/inventory/ec2.py
    $ chmod +x ec2.py
    $ export EC2_INSTANCE_FILTERS='tag:Name=FragmentorProductionDatabase'

    $ ansible -i ec2.py tag_Name_FragmentorProductionDatabase -m ping

Now you can configure the server and start and prepare the
fragmentation database. The first two plays rely on dynamic inventory
provided by the `ec2.py` script: -

    $ ansible-playbook -i ec2.py site-db-server-configure.yaml -e @parameters
    $ ansible-playbook site-db-server-configure_create-database.yaml -e @parameters

>   The Slurm compute instances may be incorrectly configured with regard to
    available memory. You need to reset the slurm manager with the correct
    memory.

The compute instances may be incorrectly configured with regard to memory.
Run the `sinfo` command to see the `MEMORY` value. If it's `1` you may need
to fix them.

    $ sinfo --exact --long -N
    [...]
    NODELIST              NODES PARTITION       STATE CPUS    S:C:T MEMORY TMP_DISK WEIGHT AVAIL_FE REASON               
    compute-dy-m4large-1      1  compute*        idle 2       2:1:1      1        0      1 dynamic, none                 
    compute-dy-m4large-2      1  compute*        idle 2       2:1:1      1        0      1 dynamic, none    

Fix the compute instance memory by running the _fix_ script.

If your compute instances have 16GiB RAM run: -

    $ ../fix-pcluster-slurm-compute-memory.sh 16

From here you should be able to run fragmentation plays, i.e. stuff like this
for a typical MolPort fragmentation run: -

    $ ansible-playbook site-standardise.yaml -e @parameters \
        -e vendor=molport \
        -e version=2020-10 

    $ ansible-playbook site-fragment.yaml -e @parameters \
        -e vendor=molport \
        -e version=2020-10

    $ ansible-playbook site-inchi.yaml -e @parameters

    $ ansible-playbook site-extract.yaml -e @parameters

---

[nextflow-pcluster]: https://github.com/InformaticsMatters/nextflow-pcluster
