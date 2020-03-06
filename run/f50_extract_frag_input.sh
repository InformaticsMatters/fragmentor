#!/bin/bash
# 
# Extract fragmentation input: 
# Purpose: Extract nonisomol data for fragmentation input
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
echo $DBHOST
echo $DATABASE
echo $PYTHONPATH/$FRAGBASEDIR/$FRAGSMIFILE
export PGPASSFILE=fragpass


echo "SMILES Extraction Starting..."
#\COPY (SELECT n.smiles FROM nonisomol n WHERE NOT EXISTS (SELECT 1 FROM edge e WHERE e.parent_id = n.id)) TO '/data/xchem/nonisomol.smi';

psql \
    -X \
    -U postgres \
    -h $DBHOST \
    --echo-all \
    --set AUTOCOMMIT=off \
    --set ON_ERROR_STOP=on \
    -c "\COPY (SELECT n.smiles FROM nonisomol n WHERE NOT EXISTS (SELECT 1 FROM edge e WHERE e.parent_id = n.id)) TO '$PYTHONPATH/$FRAGBASEDIR/$FRAGSMIFILE'" \
    $DATABASE

if [ $? -ne 0 ]; then
    echo "SMILES Extraction failed, fault:" 1>&2
    exit $?
fi

echo "SMILES Extraction Successful"
exit 0