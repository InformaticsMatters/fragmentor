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
echo $ISOINCHIFILE
echo $ISOINCHITAB
echo $INCHICHUNK
echo $VENDORPATH
source $REPPATH/$VENDORPATH/vendorparam.sh

echo "Extracting new ISOMOLs starting"
TSTART=$(date +"%T")
echo "Current time : $TSTART"

psql \
    -X \
    -U postgres \
    -h $DBHOST \
    --echo-all \
    --set AUTOCOMMIT=off \
    --set ON_ERROR_STOP=on \
    -c "\COPY (SELECT iso.smiles FROM isomol iso \
       WHERE iso.inchis IS NULL \
         AND EXISTS (SELECT 1 FROM mol_source m WHERE iso.id = m.isomol_id \
         AND m.source_id = $SOURCEID)) TO '$FRAGDATA/fragment/$ISOINCHIFILE'" \
    $DATABASE

if [ $? -ne 0 ]; then
    echo "SMILES Extraction failed, fault:" 1>&2
    exit $?
fi

echo "Calculation of Inchi starting ..."
TSTART=$(date +"%T")
echo "Current time : $TSTART"

time python -m frag.network.scripts.generate_inchi -i $FRAGPATH/fragment/$ISOINCHIFILE -o $FRAGPATH/fragment/$ISOINCHITAB -n

if [ $? -ne 0 ]; then
    echo "Fragmentation failed, fault:" 1>&2
    exit $?
fi

echo "Successfully calculated Inchi's"
TEND=$(date +"%T")
echo "Current time : $TEND"

read lines filename <<< $(wc -l $FRAGPATH/fragment/$ISOINCHITAB)
echo "lines=$lines filename=$filename"

cat $FRAGPATH/fragment/$ISOINCHITAB | split -d -l $INCHICHUNK - $FRAGPATH/fragment/isochunk_
for f in $FRAGPATH/fragment/isochunk*; do

    echo "Processing Filename $f"

    echo "Create table i_iso_inchi"

    psql \
        -X \
        -U postgres \
        -h $DBHOST \
        -f sql/f75_create_i_iso_inchi.sql \
        --echo-all \
        --set AUTOCOMMIT=on \
        --set ON_ERROR_STOP=on \
        $DATABASE

    if [ $? -ne 0 ]; then
        echo "Create table i_iso_inchi failed, fault:" 1>&2
        exit $?
    fi


    echo "Load table i_iso_inchi from inch TAB file"

    psql \
        -X \
        -U postgres \
        -h $DBHOST \
        --echo-all \
        --set AUTOCOMMIT=on \
        --set ON_ERROR_STOP=on \
        -c "\COPY i_iso_inchi(smiles, ninchik, ninchis) FROM '$f' DELIMITER E'\t' CSV;" \
        $DATABASE

    if [ $? -ne 0 ]; then
        echo "Load ISO Inchi File failed, fault:" 1>&2
        exit $?
    fi

    echo "Load Frag database from i_iso_inchi"

    psql \
        -X \
        -U postgres \
        -h $DBHOST \
        -f sql/f75_update_iso_inchi.sql \
        --echo-all \
        --set AUTOCOMMIT=off \
        --set ON_ERROR_STOP=on \
        $DATABASE

done

rm $FRAGPATH/fragment/isochunk_*

if [ $? -ne 0 ]; then
    echo "Load ISO Inchi keys failed, fault:" 1>&2
    exit $?
fi

echo "Load ISO Inchi keys successful"
TEND=$(date +"%T")
echo "Current time : $TEND"


