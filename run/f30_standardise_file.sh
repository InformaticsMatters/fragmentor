#!/bin/bash
#
# Standardise files: 
# Purpose:
#    1. Based on parameters, read input file(s) with a file template from the input directory.
#    2. For each file, call python processes to standardise a file
#    3. Concatenate standardised output into one standardised file.
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
echo $DATAPATH/data/$VENDOR
echo $STANDPATH/standardise
echo $DATAPATH/data/$VENDOR/"$STANDINPUTFILE"

export PGPASSFILE=fragpass

echo "Starting Standardisation Process .."
TSTART=$(date +"%T")
echo "Current time : $TSTART"

time nextflow run -c $REPPATH/nextflow/nextflow.config $REPPATH/nextflow/standardizer.nf -with-report $STANDPATH/standardise/standardise_report.html -with-tower \
     --script $STANDARDISER --inputs $DATAPATH/data/$VENDOR/"$STANDINPUTFILE" --out_dir $STANDPATH/standardise --chunk_size $STANDCHUNKSIZE $@

if [ $? -ne 0 ]; then
    echo "Standardisation fault, fault:" 1>&2
    exit $?
fi

echo "Standardisation Successful"
TEND=$(date +"%T")
echo "Current time : $TEND"

exit 0


