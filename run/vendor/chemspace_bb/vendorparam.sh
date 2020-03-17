#!/bin/bash
#
# Vendor Specific Parameters for DSIP
# Purpose: Environment Parameters for Fragmentation Process
#

# Filename of file to standardise - will be a txt.gz file
VENDOR='chemspace-bb'
SOURCEID=3

# Standardisation Parameters

# Module to standardise the data
STANDARDISER='frag.standardise.scripts.chemspace_bb.standardise_chemspace_bb_pricing_compounds'

# Folder that contains the input data to standardise
STANDDATADIR='run/data/chemspace_bb'

# Filename of file to standardise - will be a txt.gz file
STANDINPUTFILE='chemspace_bb.txt'

# Folder that contains the standardised data
STANDOUTPUTDIR='run/standardised-chemspace_bb'

# Chunksize for loading standardised data tab separated files
STANDCHUNK=25000
