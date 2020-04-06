#!/bin/bash
# 
# Calculate inchi key
# Purpose: Takes a list of smiles and calculates the inchi key.
# Note that f80 and f90 can be run without this step.
#
# Parameters:
#    - See file: fragparam.sh for configuration.
#    - input file : nodes.csv
#    - output file : inchis.tab for use in f85
#
# UNDER DEVELOPMENT!!!!
#
# Author | Date    | Version
# Duncan | 04/2020 | Initial Version
#

set -e
set -u

source fragparam.sh
echo $REPPATH/$FRAGBASEDIR/$FRAGNODEFILE
echo $REPPATH/$FRAGBASEDIR/$INCHITAB

export PGPASSFILE=fragpass
echo $REPPATH

echo "Calculation of Inchi starting ..."
TSTART=$(date +"%T")
echo "Current time : $TSTART"

time python -m frag.network.scripts.generate_inchi -i $REPPATH/$FRAGBASEDIR/$FRAGNODEFILE -o $REPPATH/$FRAGBASEDIR/$INCHITAB -n -s -node

if [ $? -ne 0 ]; then
    echo "Fragmentation failed, fault:" 1>&2
    exit $?
fi

echo "Successfully calculated Inchi's"
TEND=$(date +"%T")
echo "Current time : $TEND"
exit 0



