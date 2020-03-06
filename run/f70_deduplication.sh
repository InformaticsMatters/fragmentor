#!/bin/bash
# 
# De-Duplication 
# Purpose: Combine and de-duplicate edges and nodes csv files into one single file. 
# Holder step - this is only necessary for parallel streams.
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
echo $FRAGBASEDIR

