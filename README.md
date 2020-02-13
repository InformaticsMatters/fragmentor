# Fragmentation optimisation

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
