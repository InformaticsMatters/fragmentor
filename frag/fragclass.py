#!/usr/bin/env python
# Purpose
#
# Based on build_db_from_smiles.py, this module defines classes used in the
# fragmentor module.
#
# Duncan Peacock
# February 2020

import os
from multiprocessing import Process
from frag.network.models import NodeHolder, Attr
from frag.utils.network_utils import build_network

class FragProcess(Process):
    """Class FragProcess

    Purpose:
	    Manages the fragmention of a process queue element (consisting of a number of smiles).
		Returns a NodeHolder object to the results queue.

    Parameters:

    Methods:

    """

    def __init__(self, args, process_queue, results_queue) -> None:
        '''Initialises the object.'''
        super(FragProcess, self).__init__()
        self.args = args
        self.process_queue = process_queue
        self.results_queue = results_queue

    def fragment_mol(self, smiles, verbosity=0) -> object:
        """Performs the fragmentation process for a smile..

        """

        attrs = []
        attr = Attr(smiles, ["EM"])
        attrs.append(attr)

        # Build the network
        node_holder = NodeHolder(iso_flag=False)
        node_holder = build_network(attrs, node_holder, base_dir=None, verbosity=verbosity, recurse=False)

        return node_holder

    def run(self):
        """ fragment process.

        Parameters:
            A list of SMILES to process.

        Returns:
            A list of NodeHolder objects.

        """
        print('fragment process - process id: {}'.format(os.getpid()))
        while True:
            chunk_of_smiles = self.process_queue.get()
            node_list = []
            for smiles in chunk_of_smiles:
                node_holder = self.fragment_mol(smiles, verbosity=self.args.verbosity)
                node_list.append(node_holder)
            self.results_queue.put(node_list)

from threading import Thread

class FragController(Thread):
    """Class Frag Controller

    Purpose:
	   1. Manages orchestration of the Queue(s) for the FragProcesses
	   2. Manages the Cache
	   3. Manages reading of input file and writing of results??
	   4. Runs as a seperate thread - also creates a thread for filewriting?.
    Parameters:
       Input queue (Object)
       Results queue (Object)
	   input_file_name
	   reject_file_name
	   chunck_size

    Methods:

    Requires:


    """
    node_count = 0
    edge_count = 0
    reject_count = 0
    num_skipped = 0
    num_queued = 0
    num_processed = 0

    cache = set()

    def __init__(self, args, process_queue, results_queue, f_writer) -> None:
        '''Initialises the object.'''
        super(FragController, self).__init__()
        self.args = args
        self.process_queue = process_queue
        self.results_queue = results_queue
        self.f_writer = f_writer
        self.node_count = 0
        self.edge_count = 0
        self.reject_count = 0

    def get_node_count(self) -> int:
        return self.node_count

    def get_edge_count(self) -> int:
        return self.node_count

    def get_reject_count(self) -> int:
        return self.node_count

    def write_data(self, node_holder, time_ms):
        """Write_data

        Write the data from the NodeHolder to the nodes.csv and edges.csv files

        :param node_holder: Indivdual from the result queue.
        :param time_ms: The time taken to perform the fragmentation
        :param num_children: The number of child nodes
        :param num_edges: The number of child edges (some parent->child relationships have multiple edges)

        :return: the set of smiles that needs further processing.

        """
        need_further_processing = set()

        sizes = node_holder.size()  # a tuple of the number of nodes and the number of edges
        num_children = sizes[0] - 1
        num_edges = sizes[1]

        # if no edges then this is a leaf node with no children so we just write the node
        if num_edges == 0:
            if node_holder.node_list:
                node = node_holder.node_list.pop()
                smiles = node.SMILES
                if smiles not in self.cache:
                    self.f_writer.write_node(node, time_ms, 0, 0)
                    self.node_count += 1
                    self.cache.add(node.SMILES)
        else:
            # so we have edges to process
            p_node = None
            for edge in node_holder.get_edges():
                if not p_node:
                    # this happens only for the first edge
                    p_node = edge.NODES[0]
                    if p_node.SMILES in self.cache:
                        # no need to process. return immediately with need_further_processing being empty
                        return need_further_processing
                    # so it's a new node so we must write it
                    self.f_writer.write_node(p_node, time_ms, num_children, num_edges)
                    self.node_count += 1
                    self.cache.add(p_node.SMILES)
                elif edge.NODES[0] != p_node:
                    # for all other edges check that the parent is the same
                    raise ValueError("ERROR. All edges should have the same parent SMILES", p_smiles, p_node.SMILES)

                c_node = edge.NODES[1]
                p_smiles = p_node.SMILES
                c_smiles = c_node.SMILES
                if c_smiles not in self.cache:
                    need_further_processing.add(c_smiles)
                self.f_writer.write_edge(edge)
                self.edge_count += 1

        return need_further_processing

    def process_results(self, node_holder, max_frags=0, verbosity=0):
        """Process a node_holder result.

        """

        #TODO Add time_ms back in.
        time_ms = 0
        size = node_holder.size()
        # the number of children is the number of nodes minus one (the parent)
        num_children = size[0] - 1
        if 0 < max_frags < num_children:
            # TODO Reject processing needs smiles?.
            #self.f_writer.write_reject(smiles)
            self.rejects_count += 1
            return

        # print("Handling mol {0} with {1} nodes and {2} edges".format(smiles, size[0], size[1]))

        need_further_processing = self.f_writer.write_data(self, node_holder, time_ms)
        #reprocess_count = len(need_further_processing)
        # TODO Reprocessing needs chunking.
        self.process_queue.put(need_further_processing)
        self.num_queued += 1

        node_holder = None

    def read_smiles_chunk(self, standard_file, no_of_lines) -> list:
        """Read a chunk of smiles from the standard file.

        """

        smiles_to_process = []
        lines = 0

        # Change to f.readline().
        line = standard_file.readline().strip()

        while standard_file.readline() and lines < no_of_lines:
            smiles_to_process.append(line)
            line = standard_file.readline()

        return smiles_to_process

    def run(self) -> int:
        """Fragmentation Control Thread.

        Purpose:

        Fill queue - Read smiles chunk (chunk size) until max queue.

        Until no-more-smiles (file empty and no-more-results) - while queued > received
             write chunk to queue.
             loop through results
                   Add need-more-processing to queue (initially by molecule)
             If space left
                 fill with more smiles.

        Write poison pills to close down processes

        Returns:
            the number of smiles processed.
        """

        print('fragmentation control thread start')

        standard_file = open(self.args.input, 'r')
        num_skipped = self.read_smiles_chunk(self, standard_file, self.args.skip).size()
        queued_this_time = 0

        for _ in range (self.args.max_queue):
            self.process_queue.put(self.read_smiles_chunk(self, standard_file))
            self.num_queued += 1
            queued_this_time +=1

        while self.num_queued > self.num_processed:

            #Check if something has been processed
            if not self.results_queue.empty():
                node_list = self.results_queue.get()

                # Process results
                for node_holder in node_list:
                    self.process_results(self, node_holder, max_frags=self.args.max_frag, verbosity=self.args.verbosity)
                    self.num_processed += 1
                    queued_this_time +=1

            # Enough?
            if self.args.limit and self.num_processed >= self.args.limit:
                break

            # If there are not too many items queued then refill from queue
            if queued_this_time < self.args.max_queue:
                for _ in range(self.args.max_queue-queued_this_time):
                    self.process_queue.put(self.read_smiles_chunk(self, standard_file))
                    self.num_queued += 1
                queued_this_time =0

            # Check
            if self.args.report_interval > 0 and self.num_processed % self.args.report_interval == 0:
               print("Processed mol", self.num_processed)

        print('fragmentation control thread end')

        return self.num_processed

class FileWriter:
    """Class FileWriter
    Purpose:
        Contains all the output file writing for the fragmenting process.

    Methods:

    """

    def __init__(self, args, nodes_f, edges_f, rejects_f) -> None:
        '''Initialises the object.'''
        self.args = args
        self.nodes_f = nodes_f
        self.edges_f = edges_f
        self.rejects_f = rejects_f

    def write_node(self, node, time_ms, num_children, num_edges):
        """Write to Note File.

        """
        global node_count
        # print("writing node", node.SMILES)
        self.nodes_f.write(','.join([node.SMILES, str(node.HAC), str(node.RAC), str(num_children), str(num_edges), str(time_ms)]) + '\n')

    def write_edge(self, edge):
        """Write to Edge File.

        """
        global edge_count
        # print("writing edge", edge.get_label())
        self.edges_f.write(edge.as_csv() + '\n')

    def write_reject(self, smiles):
        """Write to Reject File.

        """
        self.rejects_f.write(smiles + '\n')

