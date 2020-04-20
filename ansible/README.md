# Ansible playbooks
Install requirements: -

    $ pip install -r ../requirements.txt
    $ ansible-galaxy install -r ../requirements.yaml
    
Run with something like: -

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
