#!/bin/bash
# 
# Load Standardised Data: 
# Purpose: Load Standardised Data into ISO database from results of standardisation process
# 1. Split standardised-compunds output file into chunks of chunk_size STANDCHUNK for loading into the frag database
# 2. For each chunk copy the data into i_mols and then load it into the databse
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
echo $VENDORPATH
source $REPPATH/$VENDORPATH/vendorparam.sh
echo $STANDCHUNK
echo $SOURCEID


echo $REPPATH/$STANDOUTPUTDIR

export PGPASSFILE=fragpass

head -1 $REPPATH/$STANDOUTPUTDIR/standardised-compounds.tab > standardheader
tail -n +2 $REPPATH/$STANDOUTPUTDIR/standardised-compounds.tab | split -d -l $STANDCHUNK - chunk_

for f in chunk*; do
  cat standardheader > $REPPATH/$STANDOUTPUTDIR/standardised$f.tab
  cat $f >> $REPPATH/$STANDOUTPUTDIR/standardised$f.tab

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
      echo "Copy Standardised File failed, fault:" 1>&2
      exit $?
  fi

  psql \
      -X \
      -U postgres \
      -h $DBHOST \
      -v SOURCEID=$SOURCEID \
      -f $REPPATH/$VENDORPATH/f40_load_standardised_data.sql \
      --echo-all \
      --set AUTOCOMMIT=off \
      --set ON_ERROR_STOP=on \
      $DATABASE

  if [ $? -ne 0 ]; then
      echo "Load Standardised File failed, fault:" 1>&2
      exit $?
  fi

done

rm standardheader
rm chunk*

echo "Load Standardised Results Successful"
exit 0
