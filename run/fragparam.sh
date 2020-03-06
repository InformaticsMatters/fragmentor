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

# Standardisation Parameters 

# Module to standardise the data
STANDARDISER='frag.network.scripts.standardise_xchem_compounds'

# Folder that contains the input data to standardise
STANDDATADIR='run/data'

# Filename of file to standardise - will be a txt.gz file
STANDINPUTFILE='dsip'

# Folder that contains the standardised data
STANDOUTPUTDIR='run/tmp/standardised'

# Filenames for standardised data tab separated files
STANDOUTPUTZIP='standardised-compounds.tab.gz'
STANDOUTPUTFILE='standardised-compounds.tab'

# Fragmentation Parameters 

# Module to fragment the standardised data
FRAGMENTOR='frag.network.scripts.build_db_from_smiles'

# Filename of standardised SMILES extracted from i_mols for fragmentation
FRAGSMIFILE='dsip-standardised.smi'

# Base directory for fragmentation processing - contains input and output files
FRAGBASEDIR='run/tmp'

#NB: Maybe check for these before use.
FRAGNODEFILE='nodes.csv'
FRAGEDGEFILE='edges.csv'
