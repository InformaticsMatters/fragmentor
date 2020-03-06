#!/bin/bash
# 
# Create Fragmentation Database: 
# Purpose: Creates fragmentation company specific tables
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

echo "Creating Fragmentation Database .."

psql \
    -X \
    -U postgres \
    -h $DBHOST \
    -f f20_create_frag_database.sql \
    --echo-all \
    --set AUTOCOMMIT=on \
    --set ON_ERROR_STOP=on \
    $DATABASE

if [ $? -ne 0 ]; then
    echo "Create Fragmentation Database failed, fault:" 1>&2
    exit $?
fi

echo "Create Fragmentation Database Successful"
exit 0