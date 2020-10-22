#!/usr/bin/env bash
#
# Fix the AWS ParallelCluster slurm compute memory bug.
# See https://github.com/aws/aws-parallelcluster/issues/1517
#
# Here we expect a memory size (in GiB) and we take off 2GiB
# and update the slurm config with this new value.
#
# Alan Christie
# 22 Oct 2020

ACTUAL_MEM_G=X=${1:?"Specify memory (GiB)"}

# Update slurm conf to expect memory capacity of ComputeFleet based on input
# However we lose around 800mb on average with overhead.
# Take off 2GB to be safe
REAL_MEM=$((ACTUAL_MEM_G * 1000 - 2000))
SLURM_CONF_FILE="/opt/slurm/etc/slurm.conf"
REAL_MEM_LINE="NodeName=DEFAULT RealMemory=${REAL_MEM}"
INCLUDE_CLUSTER_LINE="include slurm_parallelcluster_nodes.conf"
# Get line of INCLUDE_CLUSTER_LINE
include_cluster_line_num=$(grep -n "${INCLUDE_CLUSTER_LINE}" "${SLURM_CONF_FILE}" | cut -d':' -f1)
# Prepend with REAL_MEM_LINE
sudo sed -i "${include_cluster_line_num}i${REAL_MEM_LINE}" /opt/slurm/etc/slurm.conf
# Restart slurm database with changes to conf file
sudo systemctl restart slurmctld

# Inform the user...
echo "ComputeFleet memory set to ${REAL_MEM} MiB"
