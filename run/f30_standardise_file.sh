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
echo $PYTHONPATH
echo $VENDORPATH
source $PYTHONPATH/$VENDORPATH/vendorparam.sh

echo $STANDARDISER
echo $PYTHONPATH/$STANDDATADIR
echo $PYTHONPATH/$STANDOUTPUTDIR
echo $PYTHONPATH/$STANDDATADIR/$STANDINPUTFILE
echo $PYTHONPATH/$STANDDATADIR/$STANDOUTPUTZIP

export PGPASSFILE=fragpass

echo "Starting Standardisation Process .."

# Remove Stadardisation Output Directory
if [ -d $PYTHONPATH/$STANDOUTPUTDIR ]; then
    echo "Standardisation Output directory $PYTHONPATH/$STANDOUTPUTDIR exists"
    echo -n "Remove directory [Y/N]? "
    read choice
    if [ $choice == "Y" ] ; then
	      echo "Removing Standardisation Output Directory"
	      rm -r $PYTHONPATH/$STANDOUTPUTDIR
	  else
        echo "Standardisation step will not process if this directory exists"
        exit 1
    fi
fi

python --version
time python -m $STANDARDISER $PYTHONPATH/$STANDDATADIR $STANDINPUTFILE $PYTHONPATH/$STANDOUTPUTDIR

if [ $? -ne 0 ]; then
    echo "Standardisation failed, fault:" 1>&2
    exit $?
fi

gunzip $PYTHONPATH/$STANDOUTPUTDIR/$STANDOUTPUTZIP
#python -m frag.network.scripts.standardise_xchem_compounds ~/data dsip ~/data/standardised
#gunzip ~/data/standardised/standardised-compounds.tab.gz

if [ $? -ne 0 ]; then
    echo "Standardisation unzip, fault:" 1>&2
    exit $?
fi

echo "Standardisation Successful"
exit 0


