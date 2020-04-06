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
VENDORPATH='run/vendor/dsip'

# Fragmentation Parameters 

# Filename of standardised SMILES extracted from i_mols for fragmentation
FRAGSMIFILE='fragmentinput.smi'

# Base directory for fragmentation processing - contains input and output files
FRAGBASEDIR='run/fragment'

# Chunk size for fragmentation processing
FRAGCHUNKSIZE=15000
# FRAGCHUNKSIZE=500000
# Maximum heavy atom count for extraction/fragmentation processing
FRAGHAC=36
# Maximum frag cycles for fragmentation processing
FRAGMAXFRAGS=12

# Nodes and Edges Files produced by fragmentation - currently hardcoded into build_db_from_smiles.
FRAGNODEFILE='nodes.csv'
FRAGEDGEFILE='edges.csv'

# Chunk size for loading edges and nodes
NODECHUNK=100000
EDGECHUNK=250000

# Module to add the inchi key - creates inchi.tab from nodes.csv.
INCHITAB='inchi.tab'
INCHICHUNK=1000000

