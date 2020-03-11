#!/bin/bash
#
# Vendor Specific Parameters for DSIP
# Purpose: Environment Parameters for Fragmentation Process
#

# Filename of file to standardise - will be a txt.gz file
VENDOR='dsip'
SOURCEID=6

# Standardisation Parameters

# Module to standardise the data
STANDARDISER='frag.standardise.scripts.dsip.standardise_xchem_compounds'

# Folder that contains the input data to standardise
STANDDATADIR='run/data/dsip'

# Filename of file to standardise - will be a txt.gz file
STANDINPUTFILE='dsip.txt'

# Folder that contains the standardised data
STANDOUTPUTDIR='run/standardised-dsip'

# Filenames for standardised data tab separated files
STANDOUTPUTZIP='standardised-compounds.tab.gz'
STANDOUTPUTFILE='standardised-compounds.tab'
