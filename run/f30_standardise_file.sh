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
echo $REPPATH/$STANDDATADIR
echo $REPPATH/$STANDOUTPUTDIR
echo $REPPATH/$STANDDATADIR/$STANDINPUTFILE

export PGPASSFILE=fragpass

echo "Starting Standardisation Process .."

if [ -f $REPPATH/$STANDOUTPUTDIR/standardised-tmp.tab ] ; then
    rm $REPPATH/$STANDOUTPUTDIR/standardised-tmp.tab
fi

for f in $REPPATH/$STANDDATADIR/$STANDINPUTFILE;
do
   time nextflow run -c $REPPATH/nextflow/nextflow.config $REPPATH/nextflow/standardizer.nf -with-docker --script $STANDARDISER --input $f --out_dir $REPPATH/$STANDOUTPUTDIR --chunk_size $STANDCHUNK $@

   if [ ! -f $REPPATH/$STANDOUTPUTDIR/standardised-tmp.tab ] ; then
       cp $REPPATH/$STANDOUTPUTDIR/standardised-compounds.tab $REPPATH/$STANDOUTPUTDIR/standardised-tmp.tab
   else
       tail -n +2 $REPPATH/$STANDOUTPUTDIR/standardised-compounds.tab >> $REPPATH/$STANDOUTPUTDIR/standardised-tmp.tab
   fi

done

rm $REPPATH/$STANDOUTPUTDIR/standardised-compounds.tab
mv $REPPATH/$STANDOUTPUTDIR/standardised-tmp.tab $REPPATH/$STANDOUTPUTDIR/standardised-compounds.tab

if [ $? -ne 0 ]; then
    echo "Standardisation fault, fault:" 1>&2
    exit $?
fi

echo "Standardisation Successful"
exit 0


