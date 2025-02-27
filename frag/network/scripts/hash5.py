import argparse
import csv
import ctypes
import gzip
import os
import time
from pathlib import Path

from frag.network.scripts import patch_inchi


def run(input, path, separator=",", header=False, patch_missing_inchi=False):

    hashseed = os.getenv('PYTHONHASHSEED')
    if hashseed != '0':
        print('PYTHONHASHSEED environment variable must be set to 0 to get determinisitic hashing.')
        exit(1)

    count = 0
    collisions = 0
    num_patched_inchi = 0
    hashu=lambda word: ctypes.c_uint64(hash(word)).value

    t0 = time.time()

    with (gzip.open(input, 'rt') if input.endswith('.gz') else open(input, 'rt')) as file:
        reader = csv.reader(file, delimiter=separator)
        for row in reader:
            if count == 0 and header:
                continue
            if count % 1000000 == 0:
                print('... processed', count, collisions, num_patched_inchi)

            if patch_missing_inchi:
                ikey = row[4]
                row = patch_inchi.patch_line(row, always=True)
                if ikey != row[4]:
                    num_patched_inchi += 1
            row[5] = '"' + row[5] + '"'
            line = ','.join(row)

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
    print('Processed', count, collisions, num_patched_inchi, 'in', (t1 - t0), 'secs')


def main():
    # run like this:
    #   PYTHONHASHSEED=0 python -m frag.network.scripts.hash5 -i /work/nodes-C1.csv.gz -o /home/timbo/hashed5-C1 -p
    parser = argparse.ArgumentParser(description="Analyse molecules")

    parser.add_argument("-i", "--input", required=True, help="input csv file")
    parser.add_argument("-o", "--output", required=True, help="output dir")
    parser.add_argument("-s", "--separator", default=",", help="separator")
    parser.add_argument("-l", "--header-line", action="store_true", help="skip the first line")
    parser.add_argument("-p", "--patch-inchi", action='store_true', help="recalculate InCHI data")

    args = parser.parse_args()

    run(args.input, args.output, separator=args.separator, header=-args.header_line, patch_missing_inchi=args.patch_inchi)


if __name__ == "__main__":
    main()
