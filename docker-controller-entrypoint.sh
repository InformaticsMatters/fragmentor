#!/bin/bash

# The controller-container entrypoint.
# The environment variable FRAGMENTOR_PLAY defines
# what the container image will do.
# It is expected to be one of the plays.
# i.e. 'standardise' for the site-standardise.yaml play.

: "${FRAGMENTOR_PLAY?Need to set FRAGMENTOR_PLAY}"
echo "+> FRAGMENTOR_PLAY is ${FRAGMENTOR_PLAY}"

# You also need AWS access keys (for bucket access)
: "${AWS_ACCESS_KEY_ID?Need to set AWS_ACCESS_KEY_ID}"
: "${AWS_SECRET_ACCESS_KEY?Need to set AWS_SECRET_ACCESS_KEY}"

# A playbook parameter file '$HOME/parameters' is expected to be mapped
# into the container. Ths file contains user-specific parameters for
# the play that is to be run.

PARAMETER_FILE=${HOME}/parameters
if [ ! -f "${PARAMETER_FILE}" ]; then
    echo "+> PARAMETER_FILE (${PARAMETER_FILE}) does not exist"
    exit 1
fi

# A Nextflow configuration file '$HOME/nextflow.config' is expected to be mapped
# into the container.

NF_CONFIG_FILE=${HOME}/nextflow.config
if [ ! -f "${NF_CONFIG_FILE}" ]; then
    echo "+> NF_CONFIG_FILE (${NF_CONFIG_FILE}) does not exist"
    exit 1
fi

# All set - run the playbook...

PLAYBOOK="site-${FRAGMENTOR_PLAY}.yaml"
echo "+> Playing ${PLAYBOOK}..."
pushd ansible || exit 1
ansible-playbook "${PLAYBOOK}" -e "@${PARAMETER_FILE}"
echo "+> Played"
