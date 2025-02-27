import argparse
import time
from pathlib import Path


digits = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd', 'e', 'f']


def run(inputs, output, sections=None):

    t0 = time.time()
    pairs = []
    for i in digits:
        for j in digits:
            pairs.append(i + j)

    num_with_dups = 0
    count = 0
    non_dups = 0
    num_inputs = 0
    num_outputs = 0

    out_dir = Path(output)
    if not out_dir.exists():
        out_dir.mkdir(parents=True)

    for pair1 in pairs if not sections else sections:
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
                            for line in file:
                                tokens = line.split(',')
                                s = tokens[0]
                                smiles[s] = line
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
                        for line in smiles.values():
                            out.write(line)

        s1 = time.time()
        print('...', pair1, 'number with duplicates =', num_with_dups, 'number inputs =', num_inputs,
              'number outputs =', num_outputs, 'time =', round(s1 - s0), 'secs')
    t1 = time.time()
    print('Processing took', round(t1 - t0), 'secs')


def main():
    #  python -m frag.network.scripts.dedup5 -i ~/hashed5-C1 ~/hashed5-C2 -o ~/hashed5-out

    parser = argparse.ArgumentParser(description="collator")

    parser.add_argument("-i", "--inputs", nargs="+", help="Input dirs containing hashed data")
    parser.add_argument("-o", "--output", help="Output dir")
    parser.add_argument("-s", "--sections", nargs="*",
                        help="Top level 2 character hashes to handle (if not specified all are handled")

    args = parser.parse_args()

    run(args.inputs, args.output, sections=args.sections)


if __name__ == "__main__":
    main()