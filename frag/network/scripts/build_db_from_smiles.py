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
import time

from frag.network.models import NodeHolder, Attr
from frag.utils.network_utils import build_network

cache = set()
node_count = 0
edge_count = 0
rejects_count = 0

base_dir = None
nodes_f = None
edges_f = None
rejects_f = None

def write_node(node, time_ms, num_children, num_edges):
    global node_count
    # print("writing node", node.SMILES)
    nodes_f.write(','.join([node.SMILES, str(node.HAC), str(node.RAC), str(num_children), str(num_edges), str(time_ms)]) + '\n')
    node_count += 1
    cache.add(node.SMILES)

def write_edge(edge):
    global edge_count
    # print("writing edge", edge.get_label())
    edges_f.write(edge.as_csv() + '\n')
    edge_count += 1

def write_reject(smiles):
    global rejects_f

    if not rejects_f:
        rejects_f = open(os.path.join(base_dir, "rejects.smi"), "w")
    rejects_f.write(smiles + '\n')


def write_data(node_holder, time_ms):
    """
    Write the data from the NodeHolder to the nodes.csv and edges.csv files
    :param node_holder:
    :param time_ms: The time taken to perform the fragmentation
    :param num_children: The number of child nodes
    :param num_edges: The number of child edges (some parent->child relationships have multiple edges)
    :return:
    """
    need_further_processing = set()

    sizes = node_holder.size() # a tuple of the number of nodes and the number of edges
    num_children = sizes[0] - 1
    num_edges = sizes[1]

    # if no edges then this is a leaf node with no children so we just write the node
    if num_edges == 0:
        if node_holder.node_list:
            node = node_holder.node_list.pop()
            smiles = node.SMILES
            if smiles not in cache:
                write_node(node, time_ms, 0, 0)
    else:
        # so we have edges to process
        p_node = None
        for edge in node_holder.get_edges():
            if not p_node:
                # this happens only for the first edge
                p_node = edge.NODES[0]
                if p_node.SMILES in cache:
                    # no need to process. return immediately with need_further_processing being empty
                    return need_further_processing
                # so it's a new node so we must write it
                write_node(p_node, time_ms, num_children, num_edges)
            elif edge.NODES[0] != p_node:
                # for all other edges check that the parent is the same
                raise ValueError("ERROR. All edges should have the same parent SMILES", p_smiles, p_node.SMILES)

            c_node = edge.NODES[1]
            p_smiles = p_node.SMILES
            c_smiles = c_node.SMILES
            if c_smiles not in cache:
                need_further_processing.add(c_smiles)
            write_edge(edge)

    return need_further_processing

def fragment_and_write(smiles, max_frags=0, verbosity=0):

    global rejects_count

    t0 = time.time()
    node_holder = fragment_mol(smiles, verbosity=verbosity)
    t1 = time.time()
    time_ms = int(round((t1 - t0) * 1000))
    size = node_holder.size()
    # the number of children is the number of nodes minus one (the parent)
    num_children = size[0] - 1
    if 0 < max_frags < num_children:
        write_reject(smiles)
        rejects_count += 1
        return
    # print("Handling mol {0} with {1} nodes and {2} edges".format(smiles, size[0], size[1]))
    need_further_processing = write_data(node_holder, time_ms)
    reprocess_count = len(need_further_processing)
    for smiles in need_further_processing:
        # print("Recursing for ", child.smiles)
        fragment_and_write(smiles, max_frags=max_frags, verbosity=verbosity)
    node_holder = None

def fragment_mol(smiles, verbosity=0):

    attrs = []
    attr = Attr(smiles, ["EM"])
    attrs.append(attr)

    # Build the network
    node_holder = NodeHolder(iso_flag=False)
    node_holder = build_network(attrs, node_holder, base_dir=None, verbosity=verbosity, recurse=False)
    # Write the data out
    # print(str(node_holder.size()))
    # for node in node_holder.node_list:
    #     print(str(node))
    # for edge in node_holder.get_edges():
    #         print(str(edge))

    return node_holder

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
    nodes_f.close()
    with open(args.input, 'r') as standard_file:

        # Process the rest of the file...
        num_skipped = 0
        num_processed = 0

        for line in standard_file:

            # Do we need to skip molecules before processing?
            if num_skipped < args.skip:
                num_skipped += 1
                continue

            line = line.strip()
            if line:
                fragment_and_write(line, max_frags=args.max_frag, verbosity=args.verbosity)

            # Enough?
            num_processed += 1
            if args.report_interval > 0 and num_processed % args.report_interval == 0:
                print("Processed mol", num_processed)
            if args.limit and num_processed >= args.limit:
                break


    edges_f.close()
    if rejects_f:
        rejects_f.close()

    print("Processed {0} molecules, wrote {1} nodes and {2} edges, {3} rejects".format(num_processed, node_count, edge_count, rejects_count))


if __name__ == "__main__":
    main()
