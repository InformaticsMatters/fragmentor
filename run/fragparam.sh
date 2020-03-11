#!/bin/bash
# 
# Fragmentation Parameters 
# Purpose: Environment Parameters for Fragmentation Process
# 

# Database Parameters - Modify to correct values

DBHOST=localhost
DATABASE=mydatabase
export PGPASSFILE=fragpass

# Environment Parameters - Modify to directory of your repository.
export PYTHONPATH=/home/duncan/Documents/dev/InfoMat/fragmentor
#export PYTHONPATH=${PYTHONPATH}:/home/[your username]/Documents/.../fragmentor

# identification of run path for vendor specific processing
VENDORPATH='run/vendor/dsip'

# Fragmentation Parameters 

# Module to fragment the standardised data
FRAGMENTOR='frag.network.scripts.build_db_from_smiles'

# Filename of standardised SMILES extracted from i_mols for fragmentation
FRAGSMIFILE='fragmentinput.smi'

# Base directory for fragmentation processing - contains input and output files
FRAGBASEDIR='run/fragment'

#NB: Maybe check for these before use.
FRAGNODEFILE='nodes.csv'
FRAGEDGEFILE='edges.csv'
