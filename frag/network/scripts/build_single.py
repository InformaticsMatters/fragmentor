#!/usr/bin/env python
#
# Build nodes and edges for a single smiles
# Example usage:
# python -m frag.network.scripts.build_single --smiles 'Oc1ccc(-c2ccccc2)cc1' --id a


import argparse
import sys

from frag.network.models import NodeHolder, Attr
from frag.utils.network_utils import build_network

def process_file(file, recurse=True, no_output=False, verbosity=0):
    with open(file, 'r') as f:

        num_processed = 0

        for line in f:
            num_processed += 1
            line = line.strip()
            if line:
                print("Processing", line)
                process_smiles(line, id=str(num_processed), recurse=recurse, no_output=no_output, verbosity=0)

def process_smiles(smiles, id='1', recurse=True, no_output=False, verbosity=0):
    attrs = []
    # print("Original SMILES: " + args.smiles)
    # mol = Chem.MolFromSmiles(args.smiles)
    # if args.standardize:
    #     mol = standardize(mol)
    #     print("Standardized SMILES: " + Chem.MolToSmiles(mol))
    # smiles = Chem.CanonSmiles(Chem.MolToSmiles(mol, isomericSmiles=True))
    # print("Canonical SMILES: " + smiles)

    attr = Attr(smiles, ["EM", id])
    attrs.append(attr)
    # Build the network
    node_holder = NodeHolder(iso_flag=False)
    max_frags = 0
    # print("Recurse:",recurse)
    node_holder = build_network(attrs, node_holder,
                                max_frags, smiles, verbosity, recurse=recurse)
    # Write the data out
    if not no_output:
        for node in node_holder.get_nodes():
            print(str(node))

        for edge in node_holder.get_edges():
            print(str(edge))

        for attr in attrs:
            print(str(attr))

    print("Number of nodes: " + str(len(node_holder.get_nodes())) + " edges: " + str(len(node_holder.get_edges())))


def main():
    # Read in a SD or SMILES file - then write out into a specified directory
    parser = argparse.ArgumentParser(
        description="Convert a SMILES to nodes, edges and attributes"
    )
    parser.add_argument("--smiles")
    parser.add_argument("--file")
    parser.add_argument("--id")
    parser.add_argument("--standardize", type=bool, default=True)
    parser.add_argument("--no-recurse", action="store_true")
    parser.add_argument("--no-output", action="store_true")

    group = parser.add_mutually_exclusive_group()
    group.add_argument("-v", dest="verbosity", action="store_const", const=1)
    group.add_argument("-vv", dest="verbosity", action="store_const", const=2)

    parser.set_defaults(verbosity=0)


    args = parser.parse_args()
    if not args.no_output:
        print(args)


    # Do we have an input and base directory?
    if args.smiles:
        process_smiles(args.smiles, no_output=args.no_output)
    elif args.file:
        process_file(args.file, no_output=args.no_output)
    else:
        print('ERROR: Must specify a smiles or a file')
        sys.exit(1)



if __name__ == "__main__":
    main()
