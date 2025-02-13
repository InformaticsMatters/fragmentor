#!/usr/bin/env bash

for f in *.tgz; do
    if [ -f "$f" ]; then
        echo "decompressing $f"
        tar xfz $f
    fi
done
