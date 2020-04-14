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
echo $FRAGDATA/fragment/$FRAGSMIFILE
echo $FRAGPATH
echo $FRAGCHUNKSIZE

export PGPASSFILE=fragpass
echo $REPPATH

echo "Fragmentation Starting ..."
TSTART=$(date +"%T")
echo "Current time : $TSTART"


#time python -m frag.network.scripts.build_db_from_smiles --input /data/xchem/nonisomol.smi --base_dir /data/xchem/
time nextflow run -c $REPPATH/nextflow/nextflow.config $REPPATH/nextflow/fragmentation.nf -with-report $FRAGPATH/fragment/frag_report.html -with-tower\
    --input $FRAGDATA/fragment/$FRAGSMIFILE --out_dir $FRAGPATH/fragment --tmp_dir $FRAGPATH/fragment --chunk_size $FRAGCHUNKSIZE --max_hac $FRAGHAC --max_frag $FRAGMAXFRAGS $@


if [ $? -ne 0 ]; then
    echo "Fragmentation failed, fault:" 1>&2
    exit $?
fi

echo "Fragmentation Successful"
TEND=$(date +"%T")
echo "Current time : $TEND"

exit 0



