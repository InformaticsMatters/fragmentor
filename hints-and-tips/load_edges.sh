#!/bin/bash
# On one particular run, the fragmentation failed because the primary key
# id for the end table exceeded the size of an integer. This script was
# used to reload the edges after the database was fixed and the partially
# loaded edges were removed
#
# Purpose: Load Edges from Fragmentation step into Frag database
# Note that this uses the deduplicated edge files (*.new) that DO NOT
# ALREADY EXIST in the database.
#
# WARNING: This script is not maintained and should be used as a useful
# template only.
#
# Author | Date    | Version
# Duncan | 02/2022 | Initial Version
#

echo "Loading Edges Starting"
TSTART=$(date +"%T")
echo "Current time : $TSTART"
export PGPASSWORD=bullfinch

for f in /data/fragmentor/enamine_1/fragment/edgechunk*.new; do

    echo "Processing Filename $f"

    psql -X -U fragmentor -h 130.246.214.154 --echo-all \
         --set AUTOCOMMIT=on --set ON_ERROR_STOP=on \
         -c "\COPY edge(parent_id, child_id, label, source_id) \
             FROM '$f' DELIMITER ',' CSV;" \
        fairmolecules

#    psql -X -U fragmentor -h 130.246.214.154 --echo-all \
#         --set AUTOCOMMIT=on --set ON_ERROR_STOP=on \
#         -c "SELECT * from mol_source where id = 13;" \
#        fairmolecules

    if [ $? -ne 0 ]; then
        echo "Load Edges results failed, fault:" 1>&2
        exit $?
    fi

    TSTART=$(date +"%T")
    echo "Current time : $TSTART"

done

echo "Loading Edges Successful"
TEND=$(date +"%T")
echo "Current time : $TEND"

exit 0
