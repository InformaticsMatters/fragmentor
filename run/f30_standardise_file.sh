#!/bin/bash
#
# Standardise files: 
# Purpose: Call python processes to standardise a file
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
echo $REPPATH
echo $VENDORPATH
source $REPPATH/$VENDORPATH/vendorparam.sh

echo $STANDARDISER
echo $REPPATH/$STANDDATADIR
echo $REPPATH/$STANDOUTPUTDIR
echo $REPPATH/$STANDDATADIR/$STANDINPUTFILE
echo $REPPATH/$STANDDATADIR/$STANDOUTPUTZIP

export PGPASSFILE=fragpass

echo "Starting Standardisation Process .."

time nextflow run -c $REPPATH/nextflow/nextflow.config $REPPATH/nextflow/standardizer.nf -with-docker --script $STANDARDISER --input $REPPATH/$STANDDATADIR/$STANDINPUTFILE --out_dir $REPPATH/$STANDOUTPUTDIR

if [ $? -ne 0 ]; then
    echo "Standardisation fault, fault:" 1>&2
    exit $?
fi

echo "Standardisation Successful"
exit 0


