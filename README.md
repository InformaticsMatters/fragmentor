# Fragmentation optimisation

[![Build Status](https://travis-ci.com/InformaticsMatters/fragmentor.svg?branch=master)](https://travis-ci.com/InformaticsMatters/fragmentor)
![GitHub tag (latest SemVer)](https://img.shields.io/github/tag/informaticsmatters/fragmentor)

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


## Running the code
```
python -m frag.network.scripts.build_db_from_smiles --input data/dsip-standardised.smi --base_dir /tmp
```
