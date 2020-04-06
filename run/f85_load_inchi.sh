#!/bin/bash
# 
# Load Inchi key and update nonisomol and ismol:
# Purpose: Load Nodes from Fragmentation step into Frag database
#
# Parameters:
#    - See file: fragparam.sh and fragpass for fragmentation configuration.
# Note that this must be run after f80 -> but could be run in parallel with f90.
#
# Author | Date    | Version
# Duncan | 04/2020 | Initial Version
#


source fragparam.sh
echo $REPPATH/$FRAGBASEDIR
echo $INCHITAB
echo $INCHICHUNK

echo "Loading Inchis Starting"
TSTART=$(date +"%T")
echo "Current time : $TSTART"

read lines filename <<< $(wc -l $REPPATH/$FRAGBASEDIR/$INCHITAB)
echo "lines=$lines filename=$filename"

cat $REPPATH/$FRAGBASEDIR/$INCHITAB | split -d -l $INCHICHUNK - $REPPATH/$FRAGBASEDIR/inchchunk_
for f in $REPPATH/$FRAGBASEDIR/inchchunk*; do

    echo "Processing Filename $f"

    echo "Create table i_inchi"

    psql \
        -X \
        -U postgres \
        -h $DBHOST \
        -f f85_create_i_inchi.sql \
        --echo-all \
        --set AUTOCOMMIT=on \
        --set ON_ERROR_STOP=on \
        $DATABASE

    if [ $? -ne 0 ]; then
        echo "Create table i_inchi failed, fault:" 1>&2
        exit $?
    fi


    echo "Load table i_inchi from inch TAB file"

    psql \
        -X \
        -U postgres \
        -h $DBHOST \
        --echo-all \
        --set AUTOCOMMIT=on \
        --set ON_ERROR_STOP=on \
        -c "\COPY i_inchi(smiles, sinchik, sinchis, ninchik, ninchis) FROM '$f' DELIMITER E'\t' CSV;" \
        $DATABASE

    if [ $? -ne 0 ]; then
        echo "Load Inchi File failed, fault:" 1>&2
        exit $?
    fi

    echo "Load Frag database from i_inchi"

    psql \
        -X \
        -U postgres \
        -h $DBHOST \
        -f f85_update_inchi.sql \
        --echo-all \
        --set AUTOCOMMIT=off \
        --set ON_ERROR_STOP=on \
        $DATABASE

done

rm $REPPATH/$FRAGBASEDIR/inchchunk_*

if [ $? -ne 0 ]; then
    echo "Load Inchi keys failed, fault:" 1>&2
    exit $?
fi

echo "Load Inchi keys successful"
TEND=$(date +"%T")
echo "Current time : $TEND"

