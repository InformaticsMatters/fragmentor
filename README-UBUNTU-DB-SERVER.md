# Setting up an ubuntu (OpenStack)) database server

The demonstration server is based on the following hardware: -

- 124 cores
- 490G RAM
- 3.1 TB root volume
- ubuntu 20.04
- Optional additional volumes for `pgwal`, `pgdata` and `pgbackup`. Attached to the
  server, but mount points and mounts will be handled by an Ansible playbook.
- A suitable SSH public key installed

From this point we assume that you have access to the cluster bastion and the server.
One example is to use a bastion and jump-host via your local `~/.ssh/config`
configuration file. An example is illustrated below: -

```
Host playground-bastion
  User ubuntu
  Hostname 130.246.81.10
  IdentityFile ~/.ssh/stfc-abc-1
  ForwardAgent yes
  ServerAliveInterval 30

Host fragmentor-db
  User ubuntu
  Hostname 192.168.34.14
  IdentityFile ~/.ssh/stfc-abc-1
  ForwardAgent yes
  ServerAliveInterval 30
  ProxyJump ubuntu@playground-bastion
  ProxyCommand ssh -W %h:%p playground-bastion
```

This allows you to connect to either machine with: -

```bash
ssh playground-bastion
ssh fragmentor-db
```

## Configuring the ubuntu database server

**SSH Key**

You will need to be able to SSH to the DB host, so you might want to create an SSH key
on the bastion to allow this if you do not have one.

```bash
ssh-keygen -t ed25519 -C "fragementor@informaticsmatters.com"
```

**Docker**

The DB server will need docker. To install on ubuntu you can follow the
[Digital Ocean] instructions. Today this should also install "docker compose" (V2),
expected by our latest playbooks.

**Fragmentor User**

You will need a user on the DB server (the `db_user_account` in our playbooks)
For the `production` variables in the playbook this is `fragmentor` (by default).
So create a user and add them to the docker group so that they can run
docker commands (the user's password is not important as the playbooks
will become the user through root privilege): -

  sudo adduser fragmentor
  sudo usermod -aG docker fragmentor

**Security groups (PostgreSQL)**

You must ensure that the network Security Groups include one
that permits PostgreSQL connections from your compute cluster
(i.e. port 5432 from anywhere?), creating a new one if required.

---

[digital ocean]: https://www.digitalocean.com/community/tutorials/how-to-install-and-use-docker-on-ubuntu-20-04
