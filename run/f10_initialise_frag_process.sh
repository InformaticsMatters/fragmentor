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

echo $PYTHONPATH

# Remove Stadardisation Output Directory
if [ -d $PYTHONPATH/run ]; then
    echo "Python path set correctly"
else
    echo "Python path must be set to the base directory for the repository"
    echo "Processing cannot proceed if this is not set correctly - see fragparam.sh file in run folder"
    echo "Initialisation Failed"
    exit 1
fi

echo "Initialisation Complete"
exit 0