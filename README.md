# Fragmentation Optimisation

[![Build Status](https://travis-ci.com/InformaticsMatters/fragmentor.svg?branch=master)](https://travis-ci.com/InformaticsMatters/fragmentor)
![GitHub tag (latest SemVer)](https://img.shields.io/github/tag/informaticsmatters/fragmentor)

Optimisation of fragmentation process through the use of a postgres database to store already fragmented data. This will allow delta changes to an existing database rather than having to completely re-fragment the input files - speeding up the loading of molecules to the Neo4j database.

Currently contains: -

- Scripts to create a postgres database schema and run an end-to-end standardisation and fragmentation process to populate the database.
- Processing is run using control parameters
- Standardise and fragmentation code based on the Fragalysis Repository 
- Nextflow scripts to control cluster for Standardisation and Fragmentation steps.
- Parameter controlled chunking of input files at various stages to control throughput to the sql database.  

The new process is under development so beware of large changes. It is intended that control will ultimately be through ansible playbooks.

## Using a conda environment

Create the environment:
```
conda env create -f environment.yml
```

Activate the environment:
```
conda activate fragmentor
```

Removing the environment:
```
conda env remove --name fragmentor
```
## Creating the Postgres Database

This is based on a docker container. Precise commands to follow.

## New Process Description

The sequence diagram below shows the basic steps in the new end-to-end fragmentation process including a fragmentation database called FairMolecules. The advantage of the database approach over the current process is that each time a new dataset of molecules is provided by the vendor, the relatively lightweight standardisation step must still be performed - but only new molecules will  have to go through the fragmentation step. As this is hardware intensive, large time/cost savings should be possible. 

The new process under development is currently based on bash scripts (labelled with names f10 -> f90), but will be automated in Ansible once the major technical issues have been fully resolved. Testing has been performed for vendors dsip, chemspace and molport (still loading edges at the time of writing). The load to Neo4j is included for completeness, but has not yet been tested.  

- The process is run by an operator actor. The operator has to configure the process and place the input files in the correct location. A large proportion of this is intended to be automated.
- The Controller is currently a host machine with the fragmentor installed and from where all database-related scripts are run.
- The Cluster is a Cluster-group with a head node also with the fragmentor installed and from where nextflow-related scripts (standardiser and fragmentation) are run.  

![Fragmentor Sequence Diagram](images/FragmentorSequence.png)



## FairMolecules Database SChema

The diagram below shows the FairMolcules database schema: 

![FairMolecules Database](images/FairMoleculesDatabase.png)



Tables beginning with “i_” are used in the loading process. Notice that i_mols is vendor specific. Note that the processing to fill the inchi table is still under development.


## Creating the database

Navigate to the run directory

```
$ p10_create_frag_database.sh
```



## Running the code

Navigate to the run directory

The parameters are in files:
- run/fragparam.sh (contains general parameters)
- run/vendor/*vendor-name*/vendorparam.sh.

Specifically, the vendor path and input file name should be set.

```
$ ./f10_initialise_frag_process.sh 
$ ./f30_standardise_file.sh 
$ time ./f40_load_standardised_data.sh > .f40.log
$ ./f50_extract_frag_input.sh
$ ./f60_fragmentation.sh
$ time ./f80_load_nodes.sh > .f80.log
$ time ./f90_load_edges.sh > .f90.log
```


