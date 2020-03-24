#!/bin/bash
# 
# Fragmentation Parameters 
# Purpose: Environment Parameters for Fragmentation Process
# 

# Database Parameters - Modify to correct values

DBHOST=localhost
DATABASE=mydatabase
export PGPASSFILE=fragpass

# Environment Parameters - Repository path needs to be changed if running the bash scripts from a different directory.
export REPPATH=..

# identification of run path for vendor specific processing
# Supported vendors
#    dsip
#    chemspace_bb
#    molport
#
VENDORPATH='run/vendor/molport'

# Fragmentation Parameters 

# Module to fragment the standardised data
FRAGMENTOR='frag.network.scripts.build_db_from_smiles'

# Filename of standardised SMILES extracted from i_mols for fragmentation
FRAGSMIFILE='fragmentinput.smi'

# Base directory for fragmentation processing - contains input and output files
FRAGBASEDIR='run/fragment'

# Chunk size for fragmentation processing
FRAGCHUNKSIZE=2500
# FRAGCHUNKSIZE=500000
# Maximum heavy atom count for extraction/fragmentation processing
FRAGHAC=36
# Maximum frag cycles for fragmentation processing
FRAGMAXFRAGS=12

#NB: Maybe check for these before use.
FRAGNODEFILE='nodes.csv'
FRAGEDGEFILE='edges.csv'
