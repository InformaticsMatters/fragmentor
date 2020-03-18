"""
generate_inchi.py

Takes SMILLES as inputs and outputs the InChi and InChi Key for that molecule.
Can generate standard or non-standard InChi and the non-standard flags to use can be specified.

Example uses:

1. Generating just standard InChi (or non-standard if you instead use the -n flag).

python -m frag.network.scripts.InChi -i data/dsip-standardised.smi -o /tmp/dsip-inchi.txt -s

In this case the output is a set of lines like:
SMILES\tInChiString\tInchiKey

2. Generating standard and non-standard InChi.

python -m frag.network.scripts.generate_inchi -i data/dsip-standardised.smi -o /tmp/dsip-inchi.txt -s -n

In this case the output is a set of lines like:
SMILES\tStanadardInChiString\tStanadardInChiStringInchiKey\tNonStanadardInChiString\tNonStanadardInChiStringInchiKey

"""

import argparse, logging

from frag.utils import standardise_utils

from rdkit import Chem

logger = logging.getLogger(__name__)

def process(input, output, gen_std, flags=''):
    global logger

    print("Using flags", flags)

    with open(input, 'r') as infile:
        with open(output, 'w') as outfile:
            num_processed = 0
            for smiles in infile:
                num_processed += 1
                smiles = smiles.strip()
                if smiles:
                    mol = Chem.MolFromSmiles(smiles)
                    try:
                        if gen_std:
                            s_inchis, s_inchik = standardise_utils.gen_inchi(mol, '')
                        if flags:
                            n_inchis, n_inchik = standardise_utils.gen_inchi(mol, flags)

                        if gen_std and flags:
                            if s_inchis == n_inchis:
                                outfile.write("\t".join([smiles, s_inchis, s_inchik, '', '']) + '\n')
                            else:
                                outfile.write("\t".join([smiles, s_inchis, s_inchik, n_inchis, n_inchik]) + '\n')
                        elif gen_std:
                            outfile.write("\t".join([smiles, s_inchis, s_inchik]) + '\n')
                        elif flags:
                            outfile.write("\t".join([smiles, n_inchis, n_inchik]) + '\n')
                    except Exception as e:
                        logger.warning('gen_inchi exception for %s', smiles)
            logger.warning("Processed %s molecules", num_processed)


def main():
    global logger

    parser = argparse.ArgumentParser(description="Generate InChi and InChiKey for molecules")
    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--output")
    parser.add_argument("-s", "--standard", action='store_true', help="Generate non-standard InChi")
    parser.add_argument("-n", "--non-standard", action='store_true', help="Generate non-standard InChi")
    parser.add_argument("-f", "--flags", default='/SaveOpt /RecMet /FixedH', help="Non-standard InChi")

    args = parser.parse_args()

    if not args.non_standard and not args.standard:
        logger.warning("Must specify --standard or --non-standard or both")

    if args.non_standard:
        process(args.input, args.output, args.standard, flags=args.flags)
    else:
        process(args.input, args.output, args.standard)


if __name__ == '__main__':
    main()