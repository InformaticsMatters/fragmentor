#!/bin/env python

import gzip, sys

if len(sys.argv) > 1:
    file = sys.argv[1]
else:
    file = 'data/suppliermol-nodes.csv.gz'
print('Testing', file)

with gzip.open(file, 'rt') as input:
    count = 0
    errors = 0
    ids = set()
    for line in input:
        count += 1
        tokens = line.strip().split(',')
        fails = []
        if len(tokens) != 3:
            fails.append('Found {} tokens'.format(len(tokens)))
        else:
            if tokens[0] in ids:
                fails.append('Duplicate ID {}'.format(tokens[0]))
            else:
                ids.add(tokens[0])
            if len(tokens[1]) == 0:
                fails.append('Missing SMILES {}')
            if tokens[2] != 'Available':
                fails.append('Incorrect 3rd token {}'.format(tokens[2]))


        if fails:
            errors += 1
            print(count, fails, line.strip())

    print('Encountered', errors, 'errors')