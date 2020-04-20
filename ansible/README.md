# Ansible playbooks
Install requirements: -

    $ pip install -r ../requirements.txt
    $ ansible-galaxy install -r ../requirements.yaml
    
You will need to provide an OS_USERNAME (expected by host_vars)
for the DB server configuration. You will also need to ensure that the
user's ~/.ssh/id_rsa is set correctly so that Ansible can ssh
to the servers. If this works you should be OK to run the playbooks...

    $ ansible -m ping all

You could configure the production database server (a destructive action)
with something like: -

    $ ansible-playbook site-configure.yaml -e deployment=production
    
Then, run a playbook with something like: -

    $ ansible-playbook site-standardise.yaml
    
And pass variables in directly or via a file: -

    $ ansible-playbook site-standardise.yaml -e x=42
    $ ansible-playbook site-standardise.yaml -e @parameters

## Variable location
-   Role-specific variables that a user might change often in the corresponding
    `defaults/main.yaml`
-   Role-specific variables that a user might not change in the corresponding
    `vars/main.yaml`
-   Common variables in `group_vars/all.yaml`
