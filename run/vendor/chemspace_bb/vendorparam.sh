#!/bin/bash
#
# Vendor Specific Parameters for DSIP
# Purpose: Environment Parameters for Fragmentation Process
#

# Vendor
VENDOR='chemspace_bb'
SOURCEID=3

# Standardisation Parameters

# Module to standardise the data
STANDARDISER='frag.standardise.scripts.chemspace_bb.standardise_chemspace_bb_pricing_compounds'

# Filename of file to standardise - will be a txt.gz file
STANDINPUTFILE='chemspace_bb.txt'

# Chunksize for standardizing/loading standardised data tab separated files
STANDCHUNKSIZE=1500
STANDCHUNK=1000000
