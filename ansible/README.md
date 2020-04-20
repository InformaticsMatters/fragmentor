# Ansible playbooks
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

You can run fragmentation plays from the head node with something like: -

    $ ansible-playbook site-standardise.yaml
    
And pass variables in directly or via a file: -

    $ ansible-playbook site-standardise.yaml -e x=42
    $ ansible-playbook site-standardise.yaml -e @parameters

## Variable locations
-   Role-specific variables that a user might change often in the corresponding
    `defaults/main.yaml`
-   Role-specific variables that a user might not change in the corresponding
    `vars/main.yaml`
-   Common (global) variables can be found in `group_vars/all.yaml`
-   Host-specific variables can be found in `host_vars`
