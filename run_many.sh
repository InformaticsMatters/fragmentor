#!/bin/bash

#!/bin/bash
input=$1
while IFS= read -r line
do
  python -m frag.network.scripts.build_single --no-output --smiles "$line"
done < "$1"
