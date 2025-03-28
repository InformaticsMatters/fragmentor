# Takes a nodes.csv or edges.csv from fragmentation and "shards" the rows into a data hierarchy based on a python hash.
# Those shards can be used to deduplicate the nodes or edges (the same entry is guaranteed to be in the same shard)
# which means that the deduplication is parallelisable 256-fold (each top level 2 character directory).
# The first 5 characters (from 16) are used, the first two define a directory into which files names with all 5
# characters are placed. e.g.
#
# 00
#  | 00000
#  | 00001
#  | ...
# 01
#  | 01000
#  | ...
# ...
# ff
#  | ff000
#  | ...
#  | fffff
#
# Each file has the row from the input for that particular (first 5 character) hash.
#
# The hashing can be done on the first token or the whole line (the --use-first-token argument).
# Typically nodes should be hashed using the first token (SMILES) whilst edges should be hashed using the whole line.

import argparse
import csv
import ctypes
import gzip
import os
import time
from pathlib import Path


def run(input, path, mode, delimiter=",", header=False):

    hashseed = os.getenv('PYTHONHASHSEED')
    if hashseed != '0':
        print('PYTHONHASHSEED environment variable must be set to 0 to get determinisitic hashing.')
        exit(1)

    count = 0
    collisions = 0
    hashu=lambda word: ctypes.c_uint64(hash(word)).value

    t0 = time.time()

    with (gzip.open(input, 'rt') if input.endswith('.gz') else open(input, 'rt')) as file:
        if mode == 'isomol-molecule-edges':
            # here we just use the whole line and don't need to tokenise it
            for line in file:
                if count == 0 and header:
                    continue

                if count % 1000000 == 0:
                    print('... processed', count, collisions)
                h = hashu(line).to_bytes(8,"big").hex()

                write_data(h, path, line)
                count += 1
        else:
            # here we need the individual tokens
            reader = csv.reader(file, delimiter=delimiter)
            for row in reader:
                if count == 0 and header:
                    continue
                if count % 1000000 == 0:
                    print('... processed', count, collisions)

                if mode == 'nodes':
                    row[4] = ''
                    row[5] = ''
                elif mode == 'isomol-nodes':
                    row[1] = ''
                    row[2] = ''

                if mode == 'isomol-nodes':
                    line = row[0].strip()
                else:
                    line = ','.join(row)

                if mode == 'nodes' or mode == 'isomol-nodes':
                    h = hashu(row[0].strip()).to_bytes(8,"big").hex()
                else:
                    h = hashu(line).to_bytes(8,"big").hex()

                write_data(h, path, line + '\n')
                count += 1

    t1 = time.time()
    print('Processed', count, collisions, 'in', (t1 - t0), 'secs')


def write_data(h, path, line):
    p1 = h[0:2]
    d = Path(path) / p1

    if not d.is_dir():
        d.mkdir(parents=True)
    f = d / h[0:5]

    with open(f, "a") as out:
        out.write(line)


def main():
    # run like this:
    #   PYTHONHASHSEED=0 python -m frag.network.scripts.hash5 -i /work/nodes-C1.csv.gz -o /home/timbo/hashed5-C1 -m nodes
    #   PYTHONHASHSEED=0 python -m frag.network.scripts.hash5 -i /work/enamine_mferla/sep2024/edges.csv.gz -o /home/timbo/edges-xxx -m edges

    parser = argparse.ArgumentParser(description="Shard using a hash")

    parser.add_argument("-i", "--input", required=True, help="input csv file")
    parser.add_argument("-o", "--output", required=True, help="output dir")
    parser.add_argument("-s", "--delimiter", default=",", help="delimiter")
    parser.add_argument("-l", "--header-line", action="store_true", help="skip the first line")
    parser.add_argument("-m", "--mode", required=True,
                        choices=['nodes', 'edges', 'isomol-nodes', 'isomol-molecule-edges'],
                        help="which mode")

    args = parser.parse_args()

    run(args.input, args.output, args.mode, delimiter=args.delimiter, header=-args.header_line)


if __name__ == "__main__":
    main()
