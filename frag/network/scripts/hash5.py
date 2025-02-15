import argparse
import ctypes
import gzip
import os
import time
from pathlib import Path


def run(input, path, separator=",", token=None, header=False):

    hashseed = os.getenv('PYTHONHASHSEED')
    if hashseed != '0':
        print('PYTHONHASHSEED environment variable must be set to 0 to get determinisitic hashing.')
        exit(1)

    count = 0
    collisions = 0
    hashu=lambda word: ctypes.c_uint64(hash(word)).value

    t0 = time.time()

    with (gzip.open(input, 'rt') if input.endswith('.gz') else open(input, 'rt')) as file:
        if header:
            file.readline()
        for line in file:
            if count % 1000000 == 0:
                print('... processed', count, collisions)
            if token is None:
                data = line
            else:
                tokens = line.split(separator)
                data = tokens[token]
            h = hashu(data.strip()).to_bytes(8,"big").hex()

            p1 = h[0:2]
            d = Path(path) / p1

            if not d.is_dir():
                d.mkdir(parents=True)
            f = d / h[0:5]
            if f.is_file():
                collisions += 1
            with open(f, "a") as out:
                out.write(line)

            count += 1

    t1 = time.time()
    print('Processed', count, collisions, 'in', (t1 - t0), 'secs')


def main():
    # run like this:
    #   PYTHONHASHSEED=0 python frag/network/scripts/hash5.py -i /work/nodes-C2.csv.gz --separator ',' -o /home/timbo/hashed5-C2
    parser = argparse.ArgumentParser(description="Analyse molecules")

    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-s", "--separator", default=",")
    parser.add_argument("-t", "--token")
    parser.add_argument("-l", "--header-line", action="store_true")

    args = parser.parse_args()

    run(args.input, args.output, separator=args.separator, token=args.token, header=-args.header_line)


if __name__ == "__main__":
    main()
