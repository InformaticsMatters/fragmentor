#!/bin/bash

# Exit on error
set -e

# The controller-container entrypoint.
# The environment variable FRAGMENTOR_PLAY defines
# what the container image will do.
# It is expected to be one of the plays.
# i.e. 'standardise' for the site-standardise.yaml play.

: "${FRAGMENTOR_PLAY?Need to set FRAGMENTOR_PLAY}"
: "${FRAGMENTOR_NAMESPACE?Need to set FRAGMENTOR_NAMESPACE}"

echo "+> FRAGMENTOR_PLAY is ${FRAGMENTOR_PLAY}"
echo "+> FRAGMENTOR_NAMESPACE is ${FRAGMENTOR_NAMESPACE}"

# You also need AWS access keys (for bucket access)
: "${AWS_ACCESS_KEY_ID?Need to set AWS_ACCESS_KEY_ID}"
: "${AWS_SECRET_ACCESS_KEY?Need to set AWS_SECRET_ACCESS_KEY}"

# A playbook parameter file '$HOME/parameters.yaml' is expected to be mapped
# into the container. Ths file contains user-specific parameters for
# the play that is to be run.

PARAMETER_FILE=${HOME}/parameters.yaml
echo "+> PARAMETER_FILE is ${PARAMETER_FILE}"
if [ ! -f "${PARAMETER_FILE}" ]; then
    echo "+> PARAMETER_FILE does not exist"
    exit 1
fi

# A Nextflow configuration file '$HOME/nextflow.config' is expected to be mapped
# into the container.

NF_CONFIG_FILE=${HOME}/.nextflow/config
echo "+> NF_CONFIG_FILE is ${NF_CONFIG_FILE}"
if [ ! -f "${NF_CONFIG_FILE}" ]; then
    echo "+> NF_CONFIG_FILE does not exist"
    exit 1
fi

# A Kubernetes configuration file '$HOME/.kube/config' is expected to be mapped
# into the container, which we then copy so we can write to it...

KUBECONFIG_FILE=${HOME}/.kube/config
echo "+> KUBECONFIG_FILE is ${KUBECONFIG_FILE}"
if [ ! -f "${KUBECONFIG_FILE}" ]; then
    echo "+> KUBECONFIG_FILE does not exist"
    exit 1
fi

echo "+> KUBECONFIG is ${KUBECONFIG}"
echo "+> Copying ${KUBECONFIG_FILE} to ${KUBECONFIG}..."
cp "${KUBECONFIG_FILE}" "${KUBECONFIG}"

echo "+> kubectl version..."
kubectl version

# Now set default kubernetes namespace (i.e. the fragmentor)
echo "+> kubectl config set-context..."
kubectl config set-context --current --namespace="${FRAGMENTOR_NAMESPACE}"

# All set - run the playbook...

PLAYBOOK="site-${FRAGMENTOR_PLAY}.yaml"
echo "+> Playing ${PLAYBOOK}..."
pushd ansible || exit 1
ansible-playbook "${PLAYBOOK}" -e "@${PARAMETER_FILE}" \
  -e "ansible_python_interpreter=/usr/local/bin/python"
echo "+> Played"

KEEP_ALIVE_SECONDS=${KEEP_ALIVE_SECONDS:-0}
echo "+> KEEP_ALIVE_SECONDS is ${KEEP_ALIVE_SECONDS}"
echo "+> Sleeping (${KEEP_ALIVE_SECONDS})..."
sleep "${KEEP_ALIVE_SECONDS}"
echo "+> Slept"
