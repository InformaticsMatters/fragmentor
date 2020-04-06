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
echo $FRAGPATH/fragment/$FRAGSMIFILE
export PGPASSFILE=fragpass

echo $VENDORPATH
source $REPPATH/$VENDORPATH/vendorparam.sh

echo "SMILES Extraction Starting..."
TSTART=$(date +"%T")
echo "Current time : $TSTART"

#\COPY (SELECT n.smiles FROM nonisomol n WHERE NOT EXISTS (SELECT 1 FROM edge e WHERE e.parent_id = n.id)) TO '/data/xchem/nonisomol.smi';

psql \
    -X \
    -U postgres \
    -h $DBHOST \
    --echo-all \
    --set AUTOCOMMIT=off \
    --set ON_ERROR_STOP=on \
    -c "\COPY (SELECT n.smiles FROM nonisomol n \
         WHERE NOT EXISTS (SELECT 1 FROM edge e WHERE e.parent_id = n.id) \
         AND EXISTS (SELECT 1 FROM mol_source m WHERE n.id = m.nonisomol_id \
         AND m.source_id = $SOURCEID) \
         AND n.hac <= $FRAGHAC ) TO '$FRAGPATH/fragment/$FRAGSMIFILE'" \
    $DATABASE

if [ $? -ne 0 ]; then
    echo "SMILES Extraction failed, fault:" 1>&2
    exit $?
fi

echo "SMILES Extraction Successful"
TEND=$(date +"%T")
echo "Current time : $TEND"
exit 0