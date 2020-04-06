#!/bin/bash
#
# Vendor Specific Parameters for DSIP
# Purpose: Environment Parameters for Fragmentation Process
#

# Vendor
VENDOR='dsip'
SOURCEID=6

# Standardisation Parameters

# Module to standardise the data
STANDARDISER='frag.standardise.scripts.dsip.standardise_xchem_compounds'

# Filename of file to standardise - will be a txt.gz file
STANDINPUTFILE='dsip.txt'

# Chunksize for standardizing/loading standardised data tab separated files
STANDCHUNKSIZE=250
STANDCHUNK=1000000
