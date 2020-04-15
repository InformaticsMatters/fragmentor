#!/bin/bash
#
# Vendor Specific Parameters for Molport
# Purpose: Environment Parameters for Fragmentation Process
#

# Vendor
VENDOR='molport'
SOURCEID=7

# Standardisation Parameters

# Module to standardise the data
STANDARDISER='frag.standardise.scripts.molport.standardise_molport_compounds'

# Filename of file to standardise - will be a txt.gz file
STANDINPUTFILE="molport.*.txt"

# Chunksize for standardizing/loading standardised data tab separated files
STANDCHUNKSIZE=1500
STANDCHUNK=1000000
