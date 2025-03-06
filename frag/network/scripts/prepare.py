import argparse
import csv
import gzip
import time
from pathlib import Path

from frag.network.scripts import patch_inchi


digits = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd', 'e', 'f']


def run(input_dir, mode, sections=None, delimiter=',', generate_inchi=False):

    t0 = time.time()
    pairs = []
    for i in digits:
        for j in digits:
            pairs.append(i + j)

    count = 0
    num_inputs = 0
    num_outputs = 0
    num_patched_inchi = 0

    out_dir = Path(input_dir)
    if not out_dir.exists():
        out_dir.mkdir(parents=True)

    for pair1 in pairs if not sections else sections:
        d = Path(input_dir) / pair1
        if not d.is_dir():
            print('skipping', pair1, 'as no data found')
        else:
            print('processing', pair1)
            s0 = time.time()
            rows = []
            for pair2 in pairs:
                for pair3 in digits: # actually not a pair but a single digit
                    p1 = Path(pair1) / (pair1 + pair2 + pair3)
                    count += 1
                    p0 = Path(input_dir) / p1
                    if p0.is_file():
                        with open(p0, "rt") as file:
                            reader = csv.reader(file, delimiter=delimiter)
                            for row in reader:
                                rows.append(row)

            if len(rows) > 0:
                out_file = out_dir / (pair1 + '.txt.gz')
                print('writing', len(rows), 'rows to', out_file)
                with gzip.open(out_file, 'wt') as out:
                    for row in rows:
                        if mode == 'nodes':
                            _ik = 4
                            _is = 5
                        elif mode == 'isomol-nodes':
                            _ik = 1
                            _is = 2
                        if generate_inchi:
                            ikey = row[_ik]
                            row = patch_inchi.patch_line(row, _ik, _is, always=True)
                            if ikey != row[_ik]:
                                num_patched_inchi += 1
                        if row[_is] and row[_is][0] != '"':
                            row[_is] = '"' + row[_is] + '"'


                        out.write(','.join(row) + '\n')

            s1 = time.time()
            print('...', pair1, 'number inputs =', num_inputs,
                  'number outputs =', num_outputs, 'time =', round(s1 - s0), 'secs')
    t1 = time.time()
    print('Processing took', round(t1 - t0), 'secs')


def main():
    #  python -m frag.network.scripts.prepare -i ~/hashed5-out -p

    parser = argparse.ArgumentParser(description="collator")

    parser.add_argument("-i", "--inputs", help="Input dir containing combined data")
    parser.add_argument("-s", "--sections", nargs="*",
                        help="Top level 2 character hashes to handle (if not specified all are handled")
    parser.add_argument("-d", "--delimiter", default=",", help="file delimiter")
    parser.add_argument("-m", "--mode", required=True, choices=['nodes', 'edges', 'isomol-nodes'],
                        help="which mode")
    parser.add_argument("-p", "--patch-inchi", action="store_true", help="regenerate inchi for all entries")

    args = parser.parse_args()

    run(args.inputs, args.mode, sections=args.sections, delimiter=args.delimiter,
        generate_inchi=args.patch_inchi)


if __name__ == "__main__":
    main()