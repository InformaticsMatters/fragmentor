#!/bin/bash
# 
# Load Standardised Data: 
# Purpose: Load Standardised Data into ISO database from results of standardisation process
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
echo $VENDORPATH
source $REPPATH/$VENDORPATH/vendorparam.sh

echo $REPPATH/$STANDOUTPUTDIR/$STANDOUTPUTFILE

export PGPASSFILE=fragpass


echo "Creating Standardisation tables"

psql \
    -X \
    -U postgres \
    -h $DBHOST \
    -f $REPPATH/$VENDORPATH/f40_create_stand_database.sql \
    --echo-all \
    --set AUTOCOMMIT=on \
    --set ON_ERROR_STOP=on \
    $DATABASE

if [ $? -ne 0 ]; then
    echo "Create Standardisation tables failed, fault:" 1>&2
    exit $?
fi

echo "Load Standardised Results Starting ..."

#\COPY i_mols(osmiles, isosmiles, nonisosmiles, hac, cmpd_id) FROM '/data/xchem/standardised/standardised-compounds.tab' CSV DELIMITER E'\t' HEADER;
source $REPPATH/$VENDORPATH/f40_copy_standardised_data.sh

if [ $? -ne 0 ]; then
    echo "Load Standardised File failed, fault:" 1>&2
    exit $?
fi

psql \
    -X \
    -U postgres \
    -h $DBHOST \
    -f $REPPATH/$VENDORPATH/f40_load_standardised_data.sql \
    --echo-all \
    --set AUTOCOMMIT=off \
    --set ON_ERROR_STOP=on \
    $DATABASE

if [ $? -ne 0 ]; then
    echo "Update database with Standardised Results failed, fault:" 1>&2
    exit $?
fi

echo "Load Standardised Results Successful"
exit 0
