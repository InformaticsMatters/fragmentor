#!/usr/bin/env bash

for i in '0' '1' '2' '3' '4' '5' '6' '7' '8' '9' 'a' 'b' 'c' 'd' 'e' 'f'; do
  for j in '0' '1' '2' '3' '4' '5' '6' '7' '8' '9' 'a' 'b' 'c' 'd' 'e' 'f'; do
    echo "processing $i$j ..."
    python -m frag.network.scripts.molprops -i /work/04_xchem_dsip/hash5/edges/prepared/$i$j.txt.gz -o /work/04_xchem_dsip/hash5/edges/calculated/$i$j.txt.gz
  done
done