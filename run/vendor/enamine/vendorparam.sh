#!/bin/bash
#
# Vendor Specific Parameters for Enamine
# Purpose: Environment Parameters for Fragmentation Process
#

# Vendor
VENDOR='enamine'
SOURCEID=???

# Standardisation Parameters

# Module to standardise the data
STANDARDISER='frag.standardise.scripts.dsip.standardise_xchem_compounds'

# Filename of file to standardise - will be a txt.gz file
STANDINPUTFILE='enamine.txt'

# Chunksize for standardizing/loading standardised data tab separated files
STANDCHUNKSIZE=250
STANDCHUNK=1000000
