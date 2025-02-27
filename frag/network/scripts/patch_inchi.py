import argparse
import gzip
import sys

from frag.utils import standardise_utils
from rdkit import Chem
from rdkit import RDLogger

# Suppress basic RDKit logging...
RDLogger.logger().setLevel(RDLogger.ERROR)


def run(input, output, header=False):
    with (gzip.open(input, 'rt') if input.endswith('.gz') else open(input, 'rt')) as file:
        with open(output, 'wt') if output else sys.stdout as out:
            if header:
                file.readline()
            count = 0
            for line in file:
                count += 1
                tokens = line.split(',')
                tokens = patch_line(tokens)

                out.write(','.join(tokens))


def patch_line(tokens, always=False):

    if always or not tokens[4] or not tokens[5]:
        smiles = tokens[0]
        mol = Chem.MolFromSmiles(smiles)
        inchis = Chem.inchi.MolToInchi(mol, '/SaveOpt /RecMet /FixedH')
        inchik = Chem.inchi.InchiToInchiKey(inchis)
        tokens[4] = inchik
        tokens[5] = inchis

    return tokens


def main():

    parser = argparse.ArgumentParser(description="Generate InChi and InChiKey for molecules")
    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--output")

    args = parser.parse_args()

    run(args.input, args.output)


if __name__ == '__main__':
    main()
