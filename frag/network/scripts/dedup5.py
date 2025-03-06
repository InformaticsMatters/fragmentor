import argparse
import csv
import time
from pathlib import Path

from frag.network.scripts import patch_inchi


digits = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd', 'e', 'f']


def run(inputs, output, mode, sections=None, delimiter=',', generate_inchi=False, merge_flags=False):

    t0 = time.time()
    pairs = []
    for i in digits:
        for j in digits:
            pairs.append(i + j)

    sections_to_use = []
    if sections:
        for i in range(len(sections)):
            section = sections[i]
            if len(section) == 1:
                for d in digits:
                    sections_to_use.append(section + d)
            elif len(section) == 2:
                sections_to_use.append(section)
            else:
                raise ValueError("sections must has 1 or 2 characters")
        print('using sections', sections_to_use)

    num_with_dups = 0
    count = 0
    non_dups = 0
    num_inputs = 0
    num_outputs = 0
    num_patched_inchi = 0
    num_merged_flags = 0

    out_dir = Path(output)
    if not out_dir.exists():
        out_dir.mkdir(parents=True)

    for pair1 in pairs if not sections_to_use else sections_to_use:
        data_found = False
        for input in inputs:
            d = Path(input) / pair1
            if d.is_dir():
                data_found = True
        if not data_found:
            print('skipping', pair1, 'as no data found')
        else:
            print('processing', pair1)
            s0 = time.time()
            for pair2 in pairs:
                # print('processing', pair1, pair2)
                for pair3 in digits: # actually not a pair but a single digit
                    p1 = Path(pair1) / (pair1 + pair2 + pair3)
                    smiles = {}
                    num_lines = 0
                    count += 1

                    # if count % 10000 == 0:
                    #     print('... processed', count, non_dups, num_with_dups)
                    for input in inputs:
                        p0 = Path(input) / p1
                        if p0.is_file():
                            with open(p0, "rt") as file:
                                reader = csv.reader(file, delimiter=delimiter)
                                for row in reader:
                                    s = row[0]
                                    if s in smiles:
                                        if merge_flags:
                                            cur_flags = smiles[s][-1].split(';')
                                            new_flags = row[-1].split(';')
                                            cur_flags_len = len(cur_flags)
                                            if cur_flags != new_flags:
                                                if len(cur_flags) > len(new_flags):
                                                    bigger = cur_flags
                                                    smaller = new_flags
                                                else:
                                                    bigger = new_flags
                                                    smaller = cur_flags
                                                for flag in smaller:
                                                    if flag not in bigger:
                                                        bigger.append(flag)
                                                row[-1] = ';'.join(bigger)
                                                smiles[s] = row
                                                if len(bigger) != cur_flags_len:
                                                    num_merged_flags += 1

                                    else:
                                        smiles[s] = row

                                    num_lines += 1
                                    num_inputs += 1
                    num_outputs += len(smiles)
                    if len(smiles) == num_lines:
                        non_dups += 1
                    else:
                        num_with_dups += 1

                    if len(smiles) > 0:
                        d = out_dir / pair1
                        if not d.is_dir():
                            d.mkdir(parents=True)
                        # open the output file for writing
                        out_file = out_dir / p1
                        with open(out_file, 'wt') as out:
                            for row in smiles.values():
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
            print('...', pair1, 'number with duplicates =', num_with_dups, 'number inputs =', num_inputs,
                  'number outputs =', num_outputs, 'num merged flags =', num_merged_flags,
                  'time =', round(s1 - s0), 'secs')

    t1 = time.time()
    print('Processing took', round(t1 - t0), 'secs')


def main():
    #  python -m frag.network.scripts.dedup5 -i ~/hashed5-C1 ~/hashed5-C2 -o ~/hashed5-out

    parser = argparse.ArgumentParser(description="collator")

    parser.add_argument("-i", "--inputs", nargs="+", help="Input dirs containing hashed data")
    parser.add_argument("-o", "--output", help="Output dir")
    parser.add_argument("-s", "--sections", nargs="*",
                        help="Top level 2 character hashes to handle (if not specified all are handled")
    parser.add_argument("-d", "--delimiter", default=",", help="file delimiter")
    parser.add_argument("-m", "--mode", required=True, choices=['nodes', 'edges'], help="nodes or edges mode")
    parser.add_argument("-p", "--patch-inchi", action="store_true", help="regenerate inchi for all entries")
    parser.add_argument("-f", "--merge-flags", action="store_true", help="merge flags in last column")

    args = parser.parse_args()

    run(args.inputs, args.output, args.mode, sections=args.sections, delimiter=args.delimiter,
        generate_inchi=args.patch_inchi, merge_flags=args.merge_flags)


if __name__ == "__main__":
    main()