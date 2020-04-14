#!/bin/bash
# 
# Load Fragmentation results: 
# Purpose: Load Edges from Fragmentation step into Frag database
# 1. Split edges file into chunks of chunk_size EDGECHUNK for loading into the frag database
# 2. For each chunk copy the data into i_edge and then load it into the databse
#
# NOTE: f90_cpload_edges.sh should be used as it will be faster. This is an older method included for emergencies.
#
# Parameters:
#    - See file: fragparam.sh and fragpass for fragmentation configuration.
#
# Author | Date    | Version
# Duncan | 03/2020 | Initial Version
#


source fragparam.sh
echo $FRAGPATH/fragment
echo $FRAGEDGEFILE
source $REPPATH/$VENDORPATH/vendorparam.sh

echo "Loading Edges Starting"
TSTART=$(date +"%T")
echo "Current time : $TSTART"

read lines filename <<< $(wc -l $FRAGPATH/fragment/$FRAGEDGEFILE)
echo "lines=$lines filename=$filename"

cat $FRAGPATH/fragment/$FRAGEDGEFILE | split -d -l $EDGECHUNK - $FRAGPATH/fragment/edgechunk_
for f in $FRAGPATH/fragment/edgechunk*; do

    echo "Processing Filename $f"

    echo "Create table i_edge"

    psql \
        -X \
        -U postgres \
        -h $DBHOST \
        -f sql/f90_create_i_edges.sql \
        --echo-all \
        --set AUTOCOMMIT=on \
        --set ON_ERROR_STOP=on \
        $DATABASE

    if [ $? -ne 0 ]; then
        echo "Create table i_edge failed, fault:" 1>&2
        exit $?
    fi

    echo "Load table i_edge from edges CSV file"
    TSTART=$(date +"%T")
    echo "Current time : $TSTART"

    psql \
        -X \
        -U postgres \
        -h $DBHOST \
        --echo-all \
        --set AUTOCOMMIT=on \
        --set ON_ERROR_STOP=on \
        -c "\COPY i_edge(p_smiles, c_smiles, label) FROM '$f' DELIMITER ',' CSV;" \
        $DATABASE

    if [ $? -ne 0 ]; then
        echo "Load Edges File failed, fault:" 1>&2
        exit $?
    fi

     echo "Load Frag database from i_edge"
     TSTART=$(date +"%T")
     echo "Current time : $TSTART"

     psql \
         -X \
         -U postgres \
         -h $DBHOST \
         -v SOURCEID=$SOURCEID \
         -f sql/f90_load_frag_edges.sql \
         --echo-all \
         --set AUTOCOMMIT=off \
         --set ON_ERROR_STOP=on \
         $DATABASE

     if [ $? -ne 0 ]; then
         echo "Load Edges results failed, fault:" 1>&2
         exit $?
     fi

     TSTART=$(date +"%T")
     echo "Current time : $TSTART"

done

rm $FRAGPATH/fragment/edgechunk_*

echo "Loading Edges Successful"
TEND=$(date +"%T")
echo "Current time : $TEND"

exit 0

