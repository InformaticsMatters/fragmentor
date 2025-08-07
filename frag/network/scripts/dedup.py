"""
This is a very hacky script to deduplicate SMILES.
It is written specifically to handle Enamine output from Matteo.
Reads multiple input files and generates a single output file.
Any duplicates encountered are skipped.

Usage:
  python -m frag.network.scripts.dedup -i inputs*.cxsmiles -o dedup_output.cxsmiles

"""

import argparse


def run(inputs, output):
    smiles = set()
    num_dups = 0
    num_written = 0
    with open(output, "wt") as out:
        for input in inputs:
            print("Reading", input, "Currently have encountered", str(len(smiles)), "molecules and",
                  str(num_dups), "duplicates")
            with open(input, "rt") as f:
                for line in f:
                    tokens = line.split("\t")
                    smi = tokens[0]
                    if smi != "SMILES": # header line, which are repeated within the file!
                        if smi in smiles:
                            #print("Duplicate:", smi)
                            num_dups += 1
                        else:
                            smiles.add(smi)
                            out.write(line)
                            num_written += 1

    print(str(num_written), "molecules written", str(num_dups), "duplicates encountered")



def main():
    parser = argparse.ArgumentParser(description="collator")

    parser.add_argument("-i", "--inputs", nargs="+", help="Input files")
    parser.add_argument("-o", "--output", help="Output file")


    args = parser.parse_args()

    run(args.inputs, args.output)



if __name__ == "__main__":
    main()