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
echo $REPPATH/$FRAGBASEDIR/$FRAGSMIFILE
echo $REPPATH/$FRAGBASEDIR
echo $FRAGCHUNKSIZE

export PGPASSFILE=fragpass
echo $REPPATH

echo "Fragmentation Starting ..."

#time python -m $FRAGMENTOR --input $REPPATH/$FRAGBASEDIR/$FRAGSMIFILE --base_dir $REPPATH/$FRAGBASEDIR
#time python -m frag.network.scripts.build_db_from_smiles --input /data/xchem/nonisomol.smi --base_dir /data/xchem/
time nextflow run -c $REPPATH/nextflow/nextflow.config $REPPATH/nextflow/fragmentation.nf -with-report $REPPATH/$FRAGBASEDIR/frag_report.html\
    --input $REPPATH/$FRAGBASEDIR/$FRAGSMIFILE --out_dir $REPPATH/$FRAGBASEDIR --tmp_dir $REPPATH/$FRAGBASEDIR --chunk_size $FRAGCHUNKSIZE --max_hac $FRAGHAC --max_frag $FRAGMAXFRAGS $@



if [ $? -ne 0 ]; then
    echo "Fragmentation failed, fault:" 1>&2
    exit $?
fi

echo "Fragmentation Successful"
exit 0



