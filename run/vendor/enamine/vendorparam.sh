#!/bin/bash
#
# Vendor Specific Parameters for Enamine
# Purpose: Environment Parameters for Fragmentation Process
#

# Vendor
VENDOR='enamine'
SOURCEID=8

# Standardisation Parameters

# Module to standardise the data
STANDARDISER='frag.standardise.scripts.enamine.standardise_enamine_compounds'

# Filename of file to standardise - will be a txt.gz file
STANDINPUTFILE='enamine.txt'

# Chunksize for standardizing/loading standardised data tab separated files
STANDCHUNKSIZE=250
STANDCHUNK=1000000
