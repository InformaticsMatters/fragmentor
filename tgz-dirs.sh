#!/usr/bin/env bash

for f in *; do
    if [ -d "$f" ]; then
        echo "compressing $f"
        tar cfz "$f.tgz" $f && rm -rf $f
    fi
done
