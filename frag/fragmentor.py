#!/usr/bin/env python
# Purpose
#
# Based on build_db_from_smiles.py, this module builds the graph network nodes
# and edges from # the Informatics Matters 'standard' (uncompressed) file
# representation with fragmentation done over parallel processes.
#
# The output is these files:
# nodes.csv containing the molecules. The fields are: SMILES, HAC, RAC, NUM_CHILDREN, NUM_EDGES, TIME_MS
# edges.csv containing the edges. The fields are: PARENT_SMILES, CHILD_SMILES, LABEL
# rejects.smi containing the SMILES that were rejected because of the fragment count filter.
#
# Duncan Peacock
# February 2020


import argparse
import os
import sys
import time
import multiprocessing
from multiprocessing import Pool
# Local classes.
from frag.fragclass import FragProcess
from frag.fragclass import FragController
from frag.fragclass import FileWriter

base_dir = None
nodes_f = None
edges_f = None
rejects_f = None

def get_arguments():
    """Sets defaults and retrieves arguments from the command line for processing.
    Parameters:

    Returns:
        arguments/defaults.
    """

    parser = argparse.ArgumentParser(
        description="Convert un-compressed standard SMILES"
                    " for Astex Fragment network."
    )
    parser.add_argument("--input")
    parser.add_argument("--base_dir")
    parser.add_argument('-l', '--limit',
                        type=int, default=0,
                        help='Limit processing to the first N chunks,'
                             ' process all otherwise')
    parser.add_argument('-s', '--skip',
                        type=int, default=0,
                        help='Number of molecules to skip molecules'
                             ' in the input file')
    parser.add_argument('--max-frag',
                        type=int, default=0,
                        help='Limit processing to molecules with no more than'
                             ' this number of initial fragment (no limit if 0)')
    parser.add_argument('-r', '--report-interval', type=int, default=100, help='Reporting interval')
    parser.add_argument('-p', '--processes', type=int, default=4,
                        help='Number of parallel processes')

    # TODO HAC filter - discuss implementation - possibility of general "filter" file instead of active?.

    # TODO Strange bug if chunk set to greater than the file size - finishes but does not end master process.
    parser.add_argument('-c', '--chunk_size', type=int, default=50,
                        help='Size of chunk the SMILES will be grouped in to')
    parser.add_argument('-q', '--max_queue', type=int, default=50,
                        help='Limits how many new chunks are added to the queue')

    group = parser.add_mutually_exclusive_group()
    group.add_argument("-v", dest="verbosity", action="store_const", const=1)
    group.add_argument("-vv", dest="verbosity", action="store_const", const=2)

    parser.set_defaults(verbosity=0)
    return parser.parse_args()

def fragmentor(args):
    """ Main Function.
    1. Sets up parallel fragment worker processes using the FragWorker class
    2. Sets up fragmentation queue control thread using the FragControl class

    Parameters:
        Standardized filename
        File paths for base directory and output.
        Process parameters (max_frag_layers, max_queue_size, max_chunk_size)
        Reporting parameters

    Returns:
        processing details.

    """
    print('Fragment process - start')

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

    print ("No of Processes: ",args.processes)
    print ("Chunk size: ",args.chunk_size)
    print ("Max queue size (in terms of chunks): ", args.max_queue)

    # Create Directories for output files

    global base_dir
    global nodes_f
    global edges_f
    global rejects_f
    base_dir = args.base_dir

    if not os.path.isdir(base_dir):
        os.mkdir(base_dir)
    nodes_f = open(os.path.join(base_dir, "nodes.csv"), "w")
    edges_f = open(os.path.join(base_dir, "edges.csv"), "w")
    rejects_f = open(os.path.join(base_dir, "rejects.smi"), "w")

    f_writer = FileWriter(args, nodes_f, edges_f, rejects_f)

    num_processed = 0

    # Create Queues
    manager = multiprocessing.Manager()
    process_queue = manager.Queue()
    results_queue = manager.Queue()

    # Start Fragmentation Worker Processes
    pool = Pool(processes=args.processes)
    frag_processes = []
    for _ in range(args.processes):
        proc = FragProcess(
                args,
                process_queue,
                results_queue)
        frag_processes.append(proc)
        proc.start()

    t1 = time.time()

    print('Parallel Pool Created - Now starting fragment controller thread')

    # Start Fragmentation Control Thread
    frag = FragController(
                args,
                process_queue,
                results_queue,
                f_writer)
    frag.start()

    # Wait for thread to finish
    frag.join()

    print('Fragment controller thread ended - closing down')

    # Close processes when ended
    for proc in frag_processes:
        proc.terminate()

    # Close files
    nodes_f.close()
    edges_f.close()
    if rejects_f:
        rejects_f.close()

    print("Smiles: Processed {0} molecules. Wrote {1} nodes and {2} edges, {3} rejects"
          .format(frag.get_smiles_read(), f_writer.get_node_count(), f_writer.get_edge_count(), f_writer.get_reject_count()))
    print ("Fragmentation took:", time.time() - t1)


if __name__ == "__main__":
    args = get_arguments()
    fragmentor(args)

