# Ansible playbooks
Run with something like: -

    $ ansible-playbook site-configure.yaml
    
And pass variables in directly or via a file: -

    $ ansible-playbook site-configure.yaml -e x=42
    $ ansible-playbook site-configure.yaml -e "@parameters"

## Variable location
-   Role-specific variables that a user might change often in the corresponding
    `defaults/main.yaml`
-   Role-specific variables that a user might not change in the corresponding
    `vars/main.yaml`
-   Common variables in `group_vars/all.yaml`
