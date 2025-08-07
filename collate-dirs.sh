#!/usr/bin/env bash

for f in '0' '1' '2' '3' '4' '5' '6' '7' '8' '9' 'a' 'b' 'c' 'd' 'e' 'f'; do
  echo "processing $f? ..."
  cat $f?/* | gzip > $f.txt.gz
done

