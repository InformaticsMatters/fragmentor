#!/bin/bash
# 
# Fragmentation Parameters 
# Purpose: Environment Parameters for Fragmentation Process
# 

# identification of run path for vendor specific processing
# Supported vendors
#    dsip
#    chemspace_bb
#    molport
#
VENDORPATH='run/vendor/dsip'

# Database Parameters - Modify to correct values

DBHOST=localhost
DATABASE=mydatabase
export PGPASSFILE=fragpass

# Environment Parameters - Repository path needs to be changed if running the bash scripts from a different directory.
export REPPATH=..
# Environment Parameters - Processing path is where the data, standardise and fragment directories will be placed.
#export DATAPATH='/home/duncan/Documents/dev/InfoMat/run01'
export DATAPATH=..
#export STANDPATH='/home/duncan/Documents/dev/InfoMat/run01'
export STANDPATH='../testrun'
#export FRAGPATH='/home/duncan/Documents/dev/InfoMat/run01'
export FRAGPATH='../testrun'

# Fragmentation Parameters
# Filename of standardised SMILES extracted from i_mols for fragmentation
FRAGSMIFILE='nonisomol.smi'

# Chunk size for fragmentation processing
FRAGCHUNKSIZE=200
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

