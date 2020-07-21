# Ansible playbooks

>  Note that the ReadMe on the root directory fully explains the process.

Install requirements: -

    $ pip install -r ../requirements.txt
    $ ansible-galaxy install -r ../requirements.yaml
    
You will need to provide an `OS_USERNAME` (expected by host_vars)
in order to connect to the DB server configuration. You will also need to
ensure that the user's `~/.ssh/id_rsa` is set correctly so that Ansible can ssh
to the servers. If the following works you should be able to run the
project playbooks...

    $ ansible -m ping all

You could configure the production database server (a destructive action)
with something like: -

    $ ansible-playbook site-configure.yaml -e deployment=production

>   You only really need to run the `site-configure` play once.
    It configures the server with Docker and runs the designated database
    image and then formats the initial DB.

You can run standardisation plays from the head node with something like: -

    $ ansible-playbook site-standardise.yaml \
        -e deployment=production \
        -e vendor=xchem_dsip \
        -e version=v1 \
        -e runpath=/data/share-2/run01

A simple backup play can be used to copy the database files to the
backup volume in the DB server. It stops the database, copies the files
and then restarts the database: -

    $ ansible-playbook site-backup.yaml -e deployment=production
    
## Variable locations
-   Common (global) variables can be found in `group_vars/all.yaml`
-   Role-specific variables that a user might change often in the corresponding
    `defaults/main.yaml`
-   Role-specific variables that a user might not change in the corresponding
    `vars/main.yaml`
