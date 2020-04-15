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
echo $FRAGPATH/fragment
echo $INCHITAB
echo $INCHICHUNK

echo "Loading Inchis Starting"
TSTART=$(date +"%T")
echo "Current time : $TSTART"

read lines filename <<< $(wc -l $FRAGPATH/fragment/$INCHITAB)
echo "lines=$lines filename=$filename"

cat $FRAGPATH/fragment/$INCHITAB | split -d -l $INCHICHUNK - $FRAGPATH/fragment/inchchunk_
for f in $FRAGPATH/fragment/inchchunk*; do

    echo "Processing Filename $f"

    echo "Create table i_noniso_inchi"

    psql \
        -X \
        -U postgres \
        -h $DBHOST \
        -f sql/f85_create_i_noniso_inchi.sql \
        --echo-all \
        --set AUTOCOMMIT=on \
        --set ON_ERROR_STOP=on \
        $DATABASE

    if [ $? -ne 0 ]; then
        echo "Create table i_noniso_inchi failed, fault:" 1>&2
        exit $?
    fi


    echo "Load table i_noniso_inchi from inch TAB file"

    psql \
        -X \
        -U postgres \
        -h $DBHOST \
        --echo-all \
        --set AUTOCOMMIT=on \
        --set ON_ERROR_STOP=on \
        -c "\COPY i_noniso_inchi(smiles, sinchik, sinchis, ninchik, ninchis) FROM '$f' DELIMITER E'\t' CSV;" \
        $DATABASE

    if [ $? -ne 0 ]; then
        echo "Load Inchi File failed, fault:" 1>&2
        exit $?
    fi

    echo "Load Frag database from i_noniso_inchi"

    psql \
        -X \
        -U postgres \
        -h $DBHOST \
        -f sql/f85_update_noniso_inchi.sql \
        --echo-all \
        --set AUTOCOMMIT=off \
        --set ON_ERROR_STOP=on \
        $DATABASE

done

rm $FRAGPATH/fragment/inchchunk_*

if [ $? -ne 0 ]; then
    echo "Load Inchi keys failed, fault:" 1>&2
    exit $?
fi

echo "Load Inchi keys successful"
TEND=$(date +"%T")
echo "Current time : $TEND"


