#!/usr/bin/env python3
# coding=utf-8

"""A utility to merge/deduplicate shattered non-compressed (or compressed) nodes files.

The purpose is to allow distributed merge/deduplication process as duplicate lines are
guaranteed to co-exists in the same (smaller) output files. A sort'
can be run in a distributed fashion and the results processed to merge
nodes with the same smiles.

As the nodes.csv and isomol-nodes.csv files have a similar format, the script assumes:
1. A comparison is done on the first column (say this is the smiles)
2. If two consecutive rows have identical smiles then:
3. The other columns will be compared field by field. If they are the same, the first field is retained
4. If a field differs, then the second field is appended to the first field, separated by a semi-colon.
5. The merged row is written to the output file.

Duncan Peacock
July 2020
"""

import argparse
import gzip
import os
import csv
import logging

logger = logging.getLogger(__name__)

def processlines(csvinfile, csvoutfile):

    num_processed = 0
    matching_smiles = False
    stored_row = ['FirstRow']
    csvwriter = csv.writer(csvoutfile)

    for row in csv.reader(csvinfile):
        num_processed += 1

        if row[0]==stored_row[0]:
        #Matching smiles. Compare/merge rest of columns.
            matching_smiles = True
            for idx, val in enumerate(row):
                if val != stored_row[idx]:
                    # Non-matching fields are assumed to be lists that should be merged. This is done as a union of sets
                    # to eliminate duplicate values.
                    set_val = set(val.split(';'))
                    set_stored = set(stored_row[idx].split(';'))
                    if len(set_stored) > 0:
                        set_stored.update(set_val)
                        stored_row[idx] = ";".join(set_stored)
                    else:
                        stored_row[idx] = val
        else:
            if matching_smiles:
                matching_smiles = False
            if stored_row != ['FirstRow']:
                csvwriter.writerow(stored_row)
            stored_row = row

    if stored_row != ['FirstRow']:
        csvwriter.writerow(stored_row)

    logger.warning("Processed %s rows", num_processed)

def merge(input, output):
    """Given an input file, will merge any duplicate fields and write to the output file.

    :param input: The input file
    :param output: Output file to use
    """
    assert os.path.isfile(input)
    global logger

    if input.endswith('.gz'):
        with gzip.open(input, 'rt') as csvinfile, open(output, 'w') as csvoutfile:
            processlines(csvinfile, csvoutfile)
    else:
        with open(input, 'rt') as csvinfile, open(output, 'w') as csvoutfile:
            processlines(csvinfile, csvoutfile)

    csvoutfile.close()

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Merge records in input file with same ID")
    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--output")

    args = parser.parse_args()
    merge(args.input, args.output)
