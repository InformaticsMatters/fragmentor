#!/bin/bash
#
# Initialise steps for fragmentation process: 
# Purpose: Sets up directories/environment variables and so on.
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
echo "Initialising Fragmentation Process .."

echo $REPPATH
echo $VENDORPATH

TSTART=$(date +"%T")
echo "Current time : $TSTART"

# Check Python path.
if [ -d $REPPATH/run ]; then
    echo "Repository path set correctly"
else
    echo "Repository path must be set to the base directory for the repository"
    echo "Processing cannot proceed if this is not set correctly - see fragparam.sh file in run folder"
    echo "Initialisation Failed"
    exit 1
fi

# Check Vendor path.
if [ -d $REPPATH/$VENDORPATH ]; then
    echo "vendor path set correctly"
else
    echo "Vendor path must be set to the base directory for access to vendor specific parameters"
    echo "Processing cannot proceed if this is not set correctly - see fragparam.sh file in run folder"
    echo "Initialisation Failed"
    exit 1
fi

# Create Fragmentation base directory if it doesn't exist
if [ ! -d $REPPATH/$FRAGBASEDIR ]; then
    mkdir $REPPATH/$FRAGBASEDIR
fi    

# Set fragpass (database settings) security so that it will be picked up correctly later
# Note that the database in fragparam.sh must match the database in fragpass
if [ -f $REPPATH/run/fragpass ]; then
    chmod 0600 fragpass
    echo "Password file security changed"
else
    echo "fragpass file not found"
    echo "Initialisation Failed"
    exit 1
fi

echo
echo "Parameters for run are currently as follows:"
echo

source $REPPATH/$VENDORPATH/vendorparam.sh

read lines filename <<< $(wc -l $REPPATH/$STANDDATADIR/$STANDINPUTFILE)
echo "Vendor: $VENDOR"
echo "Molecules to Process (at least)          : $lines from filename: $filename"
echo "Note that there may be more that one file"
echo
echo "Chunksize for standardising              : $STANDCHUNKSIZE"
echo "Chunksize for uploading Standardised data: $STANDCHUNK"
echo "Chunksize for fragmentation              : $FRAGCHUNKSIZE"
echo "Chunksize for uploading nodes            : $NODECHUNK"
echo "Chunksize for uploading edges            : $EDGECHUNK"
echo "Chunksize for Inchi keys                 : $INCHICHUNK"


echo "Initialisation Complete"
TEND=$(date +"%T")
echo "Current time : $TEND"
exit 0
