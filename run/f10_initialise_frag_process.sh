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

echo "Initialisation Complete"
exit 0
