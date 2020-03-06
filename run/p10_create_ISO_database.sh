#!/bin/bash
# 
# Create ISO Database: 
# Purpose: Creates an ISO database from scratch and reinitialises information
#
# Parameters:
#    - See file: fragparam.sh and fragpass for fragmentation configuration.
#
# 
# Author | Date    | Version  
# Duncan | 03/2020 | Initial Version
#

set -e
set -u

source fragparam.sh
echo $DBHOST
echo $DATABASE

psql \
    -X \
    -U postgres \
    -h $DBHOST \
    -f p10_create_ISO_database.sql \
    --echo-all \
    --set AUTOCOMMIT=on \
    --set ON_ERROR_STOP=on \
    $DATABASE

if [ $? -ne 0 ]; then
    echo "Create ISO Database failed, fault:" 1>&2
    exit $?
fi

echo "Create ISO Database Successful"
exit 0    

