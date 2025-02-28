import argparse
import csv
import ctypes
import gzip
import os
import time
from pathlib import Path


def run(input, path, mode, delimiter=",", header=False, use_first_token=False, remove_inchi=False):

    hashseed = os.getenv('PYTHONHASHSEED')
    if hashseed != '0':
        print('PYTHONHASHSEED environment variable must be set to 0 to get determinisitic hashing.')
        exit(1)

    count = 0
    collisions = 0
    hashu=lambda word: ctypes.c_uint64(hash(word)).value

    t0 = time.time()

    with (gzip.open(input, 'rt') if input.endswith('.gz') else open(input, 'rt')) as file:
        reader = csv.reader(file, delimiter=delimiter)
        for row in reader:
            if count == 0 and header:
                continue
            if count % 1000000 == 0:
                print('... processed', count, collisions)

            if mode == 'nodes':
                if remove_inchi:
                    row[4] = ''
                    row[5] = ''
                else:
                    row[5] = '"' + row[5] + '"'

            line = ','.join(row)

            if use_first_token:
                h = hashu(row[0].strip()).to_bytes(8,"big").hex()
            else:
                h = hashu(line).to_bytes(8,"big").hex()

            p1 = h[0:2]
            d = Path(path) / p1

            if not d.is_dir():
                d.mkdir(parents=True)
            f = d / h[0:5]
            if f.is_file():
                collisions += 1
            with open(f, "a") as out:
                out.write(line + '\n')

            count += 1

    t1 = time.time()
    print('Processed', count, collisions, 'in', (t1 - t0), 'secs')


def main():
    # run like this:
    #   PYTHONHASHSEED=0 python -m frag.network.scripts.hash5 -i /work/nodes-C1.csv.gz -o /home/timbo/hashed5-C1 -m nodes
    #   PYTHONHASHSEED=0 python -m frag.network.scripts.hash5 -i /work/enamine_mferla/sep2024/edges.csv.gz -o /home/timbo/edges-xxx -m edges
    parser = argparse.ArgumentParser(description="Analyse molecules")

    parser.add_argument("-i", "--input", required=True, help="input csv file")
    parser.add_argument("-o", "--output", required=True, help="output dir")
    parser.add_argument("-s", "--delimiter", default=",", help="delimiter")
    parser.add_argument("-t", "--use-first-token", action="store_true", help="use first token for hashing (if not specified then whole line)")
    parser.add_argument("-l", "--header-line", action="store_true", help="skip the first line")
    parser.add_argument("-m", "--mode", required=True, choices=['nodes', 'edges'], help="nodes or edges mode")
    parser.add_argument("-r", "--remove-inchi", action="store_true", help="remove inchi columns")

    args = parser.parse_args()

    run(args.input, args.output, args.mode, delimiter=args.delimiter, header=-args.header_line,
        use_first_token=args.use_first_token, remove_inchi=args.remove_inchi)


if __name__ == "__main__":
    main()
