#!/bin/bash
# 
# Fragmentation 
# Purpose: Fragmentation of molecules
#
# Parameters:
#    - See file: fragparam.sh and fragpass for fragmentation configuration.
#
# Author | Date    | Version
# Duncan | 03/2020 | Initial Version
#

set -e
set -u

source fragparam.sh
echo $FRAGMENTOR
echo $PYTHONPATH/$FRAGBASEDIR/$FRAGSMIFILE
echo $PYTHONPATH/$FRAGBASEDIR

export PGPASSFILE=fragpass
echo $PYTHONPATH

python --version

echo "Fragmentation Starting ..."

#time python -m $FRAGMENTOR --input $PYTHONPATH/$FRAGBASEDIR/$FRAGSMIFILE --base_dir $PYTHONPATH/$FRAGBASEDIR
#time python -m frag.network.scripts.build_db_from_smiles --input /data/xchem/nonisomol.smi --base_dir /data/xchem/
nextflow run -c nextflow/nextflow.config nextflow/fragmentation.nf -with-docker --input $PYTHONPATH/$FRAGBASEDIR/$FRAGSMIFILE --out_dir $PYTHONPATH/$FRAGBASEDIR


if [ $? -ne 0 ]; then
    echo "Fragmentation failed, fault:" 1>&2
    exit $?
fi

echo "Fragmentation Successful"
exit 0



