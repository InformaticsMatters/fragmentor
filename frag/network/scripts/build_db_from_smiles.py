#!/usr/bin/env python
#
# Based on build_db.py, this module builds the graph network nodes and edges from
# the Informatics Matters 'standard' (uncompressed) file representation.
# The output is these files:
# nodes.csv containing the molecules. The fields are: SMILES, HAC, RAC, NUM_CHILDREN, NUM_EDGES, TIME_MS
# edges.csv containing the edges. The fields are: PARENT_SMILES, CHILD_SMILES, LABEL
# rejects.smi containing the SMILES that were rejected because of the fragment count filter.
#
# Tim Dudgeon
# February 2020

import argparse
import os
import sys
import collections

from frag.network.models import NodeHolder, Attr
from frag.utils.network_utils import build_network

from rdkit import Chem

class FragData():

    def __init__(self, nodes, input_smiles):
        self.parent_data = collections.OrderedDict()
        self.nodes_map = {}
        for node in nodes:
            self.nodes_map[node.SMILES] = node
        for smiles in input_smiles:
            node = self.nodes_map[smiles]
            self.parent_data[smiles] = ParentData(smiles, node.HAC, node.RAC)

    def add_edge(self, edge):
        p_smiles = edge.NODES[0].SMILES
        if p_smiles in self.parent_data:
            p_data = self.parent_data[p_smiles]
        else:
            # it's new to this group of edges
            node = self.nodes_map[p_smiles]
            p_data = ParentData(p_smiles, node.HAC, node.RAC)
            self.parent_data[p_smiles] = p_data

        p_data.add_edge(edge)

    def get_parent_data_list(self):
        """
        Get the ParentData items as a list.
        Unlike the FragData instance that list is a lightweight object suitable for pickling.
        :return:
        """
        return list(self.parent_data.values())


class ParentData():
    def __init__(self, smiles, hac, rac):
        """
        :param smiles:
        :param hac:
        :param rac:
        """
        self.smiles = smiles
        self.hac = hac
        self.rac = rac
        # children are a dict keyed by the child SMILES whose values are a list of labels
        self.children = collections.OrderedDict()
        self.edge_count = 0

    def add_edge(self, edge):
        c_smiles = edge.NODES[1].SMILES
        label = edge.get_label()
        self.edge_count += 1
        if c_smiles in self.children:
            self.children[c_smiles].append(label)
        else:
            self.children[c_smiles] = [ label ]


cache = set()
node_count = 0
edge_count = 0
rejects_count = 0

base_dir = None
nodes_f = None
edges_f = None
rejects_f = None

def write_node(smiles, hac, rac, num_children, num_edges):
    global node_count
    # print("writing node", smiles)
    # nodes_f.write(smiles + '\n')
    nodes_f.write(','.join([smiles, str(hac), str(rac), str(num_children), str(num_edges)]) + '\n')
    node_count += 1
    cache.add(smiles)

def write_edge(p_smiles, c_smiles, label):
    global edge_count
    # print("writing edge", label)
    edges_f.write(','.join([p_smiles, c_smiles, label]) + '\n')
    edge_count += 1

def write_reject(smiles):
    global rejects_f

    if not rejects_f:
        rejects_f = open(os.path.join(base_dir, "rejects.smi"), "w")
    rejects_f.write(smiles + '\n')

def write_data(parent_data):
    write_nodes(parent_data)
    return write_edges(parent_data)

def write_edges(parent_data):
    """
    Write the data to the edges.csv file
    :param parent_data: A ParentData instance

    :return:
    """
    need_further_processing = set()
    p_smiles = parent_data.smiles

    for c_smiles in parent_data.children:
        labels = parent_data.children[c_smiles]
        for label in labels:
            write_edge(p_smiles, c_smiles, label)
        if c_smiles not in cache:
            need_further_processing.add(c_smiles)

    return need_further_processing

def write_nodes(parent_data):
    """
    Write the data to the nodes.csv file
    :param parent_data:
    :return:
    """
    p_smiles = parent_data.smiles
    num_children = len(parent_data.children)
    num_edges = parent_data.edge_count
    write_node(p_smiles, parent_data.hac, parent_data.rac, num_children, num_edges)


def fragment_and_write(input_smiles, max_frags=0, max_hac=0, verbosity=0):

    #print("Fragmenting", len(input_smiles), ','.join(sorted(input_smiles)))

    global rejects_count

    if max_hac > 0:
        # note this is inefficient as we need to create the RDKit mol.
        # Better to filter the inputs.
        new_list = []
        for smiles in input_smiles:
            mol = Chem.MolFromSmiles(smiles)
            if mol.GetNumHeavyAtoms() > max_hac:
                write_reject(smiles)
                rejects_count += 1
            else:
                new_list.append(smiles)
        if not new_list:
            return
        input_smiles = new_list

    frag_data = fragment_mols(input_smiles, verbosity=verbosity)

    all_needing_further_processing = set()
    parent_smiles = set()
    data = frag_data.get_parent_data_list()

    for parent_data in data:
        parent_smiles.add(parent_data.smiles)
        if parent_data.smiles in cache:
            data.remove(parent_data)
            print("Item", parent_data.smiles, "already processed")
        else:
            num_children = len(parent_data.children)
            if 0 < max_frags < num_children:
                write_reject(parent_data.smiles)
                rejects_count += 1
                data.remove(parent_data)

    for parent_data in data:
        # print("Handling mol {0} with {1} nodes and {2} edges".format(smiles, size[0], size[1]))
        need_further_processing = write_data(parent_data)
        for item in need_further_processing:
            if item not in parent_smiles:
                all_needing_further_processing.add(item)

    for smiles in all_needing_further_processing:
        # print("Recursing for", smiles)
        # we need another check in the cache as an earlier iteration of this for loop might have handled the mol
        if smiles not in cache:
            # set max_frags and max_hac to zero as we only want filtering at the outer level.
            fragment_and_write([smiles], max_frags=0, max_hac=0, verbosity=verbosity)

    # if all_needing_further_processing:
    #     fragment_and_write(all_needing_further_processing, max_frags=0, max_hac=0, verbosity=verbosity)

def fragment_mols(input_smiles, verbosity=0):

    #print("Fragmenting", smiles)

    attrs = []
    for smiles in input_smiles:
        attr = Attr(smiles, ["EM"])
        attrs.append(attr)

    # Build the network
    # print("Processing", len(input_smiles), "mols")
    # print('INPUT ', ','.join(sorted(input_smiles)))
    node_holder = NodeHolder(iso_flag=False)
    node_holder = build_network(attrs, node_holder, base_dir=None, verbosity=verbosity, recurse=False)
    # output_smiles = [n.SMILES for n in node_holder.get_nodes()]
    # print('OUTPUT', ','.join(sorted(output_smiles)))

    frag_data = group_data(node_holder, input_smiles)
    # print("Groups:", len(frag_data.parent_data))

    return frag_data


def group_data(node_holder, input_smiles):
    """
    group the edges by parent and child SMILES

    :param node_holder:
    :param input_smiles: List of smiles that generated the node_holder
    :return: A FragData instance
    """

    frag_data = FragData(node_holder.get_nodes(), input_smiles)

    for edge in node_holder.get_edges():
        frag_data.add_edge(edge)

    return frag_data


def main():
    """Read in a 'standard' file - then write out into a specified directory
    """
    parser = argparse.ArgumentParser(
        description="Convert un-compressed standard SMILES"
                    " for Astex Fragment network."
    )
    parser.add_argument("--input")
    parser.add_argument("--base_dir")
    parser.add_argument('-l', '--limit',
                        type=int, default=0,
                        help='Limit processing to the first N molecules,'
                             ' process all otherwise')
    parser.add_argument('-s', '--skip',
                        type=int, default=0,
                        help='Number of molecules to skip molecules'
                             ' in the input file')
    parser.add_argument('--max-frag',
                        type=int, default=0,
                        help='Limit processing to molecules with no more than'
                             ' this number of initial fragment (no limit if 0)')
    parser.add_argument('--max-hac',
                        type=int, default=0,
                        help='Limit processing to molecules with no more than this heavy atoms (no limit if 0)')
    parser.add_argument('--chunk', type=int, default=1, help='Process this many smiles at a time')
    parser.add_argument('-r', '--report-interval', type=int, default=1000, help='Reporting interval')


    group = parser.add_mutually_exclusive_group()
    group.add_argument("-v", dest="verbosity", action="store_const", const=1)
    group.add_argument("-vv", dest="verbosity", action="store_const", const=2)

    parser.set_defaults(verbosity=0)
    args = parser.parse_args()

    # Do we have an input and base directory?
    if not args.input:
        print('ERROR: Must specify an input')
        sys.exit(1)
    if not os.path.isfile(args.input):
        print('ERROR: input (%s) does not exist' % args.input)
        sys.exit(2)
    if not args.base_dir:
        print('ERROR: Must specify a base directory')
        sys.exit(3)

    global base_dir
    global nodes_f
    global edges_f
    global rejects_f
    base_dir = args.base_dir

    if not os.path.isdir(base_dir):
        os.mkdir(base_dir)

    nodes_f = open(os.path.join(base_dir, "nodes.csv"), "w")
    edges_f = open(os.path.join(base_dir, "edges.csv"), "w")

    with open(args.input, 'r') as standard_file:

        # Process the rest of the file...
        num_skipped = 0
        num_processed = 0
        chunk = []

        for line in standard_file:

            num_processed += 1
            # Do we need to skip molecules before processing?
            if num_skipped < args.skip:
                num_skipped += 1
                continue

            line = line.strip()
            if line and line not in cache:
                chunk.append(line)
                if num_processed % args.chunk == 0:
                    # print("Sending chunk of", len(chunk))
                    fragment_and_write(chunk, max_frags=args.max_frag, max_hac=args.max_hac, verbosity=args.verbosity)
                    chunk = []

            # Enough?
            if args.report_interval > 0 and num_processed % args.report_interval == 0:
                print("Processed mol", num_processed)
            if args.limit and num_processed >= args.limit:
                if chunk:
                    fragment_and_write(chunk, max_frags=args.max_frag, max_hac=args.max_hac, verbosity=args.verbosity)
                print("Terminating after mol", num_processed)
                break

        if chunk:
            fragment_and_write(chunk, max_frags=args.max_frag, max_hac=args.max_hac, verbosity=args.verbosity)

    nodes_f.close()
    edges_f.close()
    if rejects_f:
        rejects_f.close()

    print("Processed {0} molecules, wrote {1} nodes and {2} edges, {3} rejects".format(num_processed, node_count, edge_count, rejects_count))


if __name__ == "__main__":
    main()
