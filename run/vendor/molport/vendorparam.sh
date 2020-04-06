#!/bin/bash
#
# Vendor Specific Parameters for Molport
# Purpose: Environment Parameters for Fragmentation Process
#

# Filename of file to standardise - will be a txt.gz file
VENDOR='molport'
SOURCEID=7

# Standardisation Parameters

# Module to standardise the data
STANDARDISER='frag.standardise.scripts.molport.standardise_molport_compounds'

# Folder that contains the input data to standardise
STANDDATADIR='run/data/molport'

# Filename of file to standardise - will be a txt.gz file
STANDINPUTFILE="molport.*.txt"

# Folder that contains the standardised data
STANDOUTPUTDIR='run/standardised-molport'

# Chunksize for standardizing/loading standardised data tab separated files
STANDCHUNKSIZE=1500
STANDCHUNK=1000000
