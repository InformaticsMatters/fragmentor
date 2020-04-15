#!/bin/bash
# 
# Load Fragmentation results: 
# Purpose: Load Nodes from Fragmentation step into Frag database
#
# Parameters:
#    - See file: fragparam.sh and fragpass for fragmentation configuration.
#
# Author | Date    | Version
# Duncan | 03/2020 | Initial Version
#


source fragparam.sh
echo $REPPATH/$FRAGBASEDIR
echo $FRAGNODEFILE
source $REPPATH/$VENDORPATH/vendorparam.sh

echo "Loading Nodes Starting"
TSTART=$(date +"%T")
echo "Current time : $TSTART"

read lines filename <<< $(wc -l $FRAGPATH/fragment/$FRAGNODEFILE)
echo "lines=$lines filename=$filename"

cat $FRAGPATH/fragment/$FRAGNODEFILE | split -d -l $NODECHUNK - $FRAGPATH/fragment/nodechunk_
for f in $FRAGPATH/fragment/nodechunk*; do

    echo "Processing Filename $f"

    echo "Create table i_node"

    psql \
        -X \
        -U postgres \
        -h $DBHOST \
        -f sql/f80_create_i_nodes.sql \
        --echo-all \
        --set AUTOCOMMIT=on \
        --set ON_ERROR_STOP=on \
        $DATABASE

    if [ $? -ne 0 ]; then
        echo "Create table i_node failed, fault:" 1>&2
        exit $?
    fi


    echo "Load table i_node from nodes CSV file"

    psql \
        -X \
        -U postgres \
        -h $DBHOST \
        --echo-all \
        --set AUTOCOMMIT=on \
        --set ON_ERROR_STOP=on \
        -c "\COPY i_node(smiles, hac, rac, ring_smiles, child_count, edge_count) FROM '$f' DELIMITER ',' CSV;" \
        $DATABASE

    if [ $? -ne 0 ]; then
        echo "Load Nodes File failed, fault:" 1>&2
        exit $?
    fi

    echo "Load Frag database from i_nodes"

    psql \
        -X \
        -U postgres \
        -h $DBHOST \
        -v SOURCEID=$SOURCEID \
        -f sql/f80_load_frag_nodes.sql \
        --echo-all \
        --set AUTOCOMMIT=off \
        --set ON_ERROR_STOP=on \
        $DATABASE

done

rm $FRAGPATH/fragment/nodechunk_*

if [ $? -ne 0 ]; then
    echo "Load Nodes results failed, fault:" 1>&2
    exit $?
fi

echo "Load nodes successful"
TEND=$(date +"%T")
echo "Current time : $TEND"


