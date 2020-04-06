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

3. Adding standard InChi to a nodes.csv file generated from build_db_from_smiles.

python -m frag.network.scripts.generate_inchi -i data/dsip-standardised.smi -o /tmp/dsip-inchi.txt -snode -n

In this case the output will be a TAB (changed from csv) file with a set of lines like:
SMILES\thac\trac\tchild_count\tedge_count\tInChiString\tInchiKey
"""

import argparse
import logging
import csv


from frag.utils import standardise_utils
from rdkit import Chem
from rdkit import RDLogger

logger = logging.getLogger(__name__)

# Suppress basic RDKit logging...
RDLogger.logger().setLevel(RDLogger.ERROR)

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

def process_node(input, output, gen_std, flags=''):
    """Reads a nodes.csv file, adds the inchi key and writes as nodes.tab.
    Parameters:
        nodes.csv

    Returns:
        inchi.tab.
    """
    global logger

    with open(input, 'r') as csvinfile, open(output, 'w') as taboutfile:
        num_processed = 0
        errors = 0

        fieldnames = ['smiles', 'hac','rac','child_count','edge_count']
        reader = csv.DictReader(csvinfile,fieldnames=fieldnames[0:5])

        for row in reader:
            num_processed += 1
            s_inchis=''
            s_inchik=''
            smiles = row['smiles'].strip()
            if smiles:
                mol = Chem.MolFromSmiles(smiles)
                try:
                    if gen_std:
                        s_inchis, s_inchik = standardise_utils.gen_inchi(mol, '')
                    if flags:
                        n_inchis, n_inchik = standardise_utils.gen_inchi(mol, flags)

                    if gen_std and flags:
                        if s_inchis == n_inchis:
                            taboutfile.write("\t".join([smiles, s_inchis, s_inchik, '', '']) + '\n')
                        else:
                            taboutfile.write("\t".join([smiles, s_inchis, s_inchik, n_inchis, n_inchik]) + '\n')
                    elif gen_std:
                        taboutfile.write("\t".join([smiles, s_inchis, s_inchik]) + '\n')
                    elif flags:
                        taboutfile.write("\t".join([smiles, n_inchis, n_inchik]) + '\n')
                except Exception as e:
                    errors += 1
                    logger.warning('gen_inchi exception for %s', smiles)
        logger.warning("Processed %s molecules", num_processed)
        logger.warning("Errors Generating %s Keys", errors)

def main():
    global logger

    parser = argparse.ArgumentParser(description="Generate InChi and InChiKey for molecules")
    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--output")
    parser.add_argument("-s", "--standard", action='store_true', help="Generate standard InChi")
    parser.add_argument("-n", "--non-standard", action='store_true', help="Generate non-standard InChi")
    parser.add_argument("-f", "--flags", default='/SaveOpt /RecMet /FixedH', help="Non-standard InChi")
    parser.add_argument("-node", "--node-file", action='store_true', help="Input file is a nodes.csv file")

    args = parser.parse_args()

    if not args.non_standard and not args.standard:
        logger.warning("Must specify --standard or --non-standard")

    if args.node_file:
        process_node(args.input, args.output, args.standard, flags=args.flags)
    elif args.non_standard:
        process(args.input, args.output, args.standard, flags=args.flags)
    else:
        process(args.input, args.output)

if __name__ == '__main__':
    main()