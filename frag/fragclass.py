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

    """

    def __init__(self, args, process_queue, results_queue) -> None:
        '''Initialises the object.'''
        super(FragProcess, self).__init__()
        self.args = args
        self.process_queue = process_queue
        self.results_queue = results_queue

    def fragment_mol(self, smiles, verbosity=0) -> object:
        """Performs the fragmentation process for a SMILES.

        Returns:
           NodeHolder object.

        """

        attrs = []
        attr = Attr(smiles, ["EM"])
        attrs.append(attr)
        #print('fragment smiles: {}'.format(smiles))
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
            if not self.process_queue.empty():
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

    """

    #Queue control
    num_queued = 0
    num_processed = 0
    num_requeued = 0
    queued_this_time = 0

    smiles_read = 0
    cache = set()

    def __init__(self, args, process_queue, results_queue, f_writer) -> None:
        '''Initialises the object.'''
        super(FragController, self).__init__()
        self.args = args
        self.process_queue = process_queue
        self.results_queue = results_queue
        self.f_writer = f_writer

    def get_num_processed(self) -> int:
        return self.num_processed

    def get_smiles_read(self) -> int:
        return self.smiles_read

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

        return need_further_processing

    def process_results(self, node_holder, max_frags=0, verbosity=0):
        """Processes a NodeHolder object from the results queue.

        """

        #TODO Add time_ms back in - has to be returned with node holders.
        time_ms = 0
        size = node_holder.size()
        # the number of children is the number of nodes minus one (the parent)
        num_children = size[0] - 1
        if 0 < max_frags < num_children:
            # TODO Reject processing needs smiles?.
            # self.f_writer.write_reject(smiles)
            return

        # print("Handling mol {0} with {1} nodes and {2} edges".format(smiles, size[0], size[1]))

        need_further_processing = self.write_data(node_holder, time_ms)
        #reprocess_count = len(need_further_processing)
        # TODO Reprocessing needs chunking.
        if len(need_further_processing) > 0:
           self.process_queue.put(need_further_processing)
           self.num_requeued += 1
           self.queued_this_time += 1

        node_holder = None

    def read_smiles_chunk(self, standard_file, no_of_lines) -> list:
        """Read a chunk of SMILES from the standard file.

        Returns:
            A list of SMILES to add to the process queue.

        """
        smiles_to_process = []

        while True and len(smiles_to_process) < no_of_lines:
            # read a single line
            line = standard_file.readline().strip()
            if line:
                smiles_to_process.append(line)
            else:
                break

        self.smiles_read += len(smiles_to_process)
        print ("Smiles read :", self.smiles_read)
        return smiles_to_process

    def run(self) -> int:
        """Fragmentation Control Thread.

        Fill queue - Read smiles chunk (chunk size) until max queue.
        while queued > received
             loop through results
                 Add Nodes that need-more-processing to queue
             If space left
                 Fill queue - Read smiles chunk (chunk size) until max queue.
        """

        print('Fragmentation Control Thread Start')

        standard_file = open(self.args.input, 'r')
        smiles_to_process = True

        if self.args.skip > 0:
            num_skipped = self.read_smiles_chunk(standard_file, self.args.skip).size()
        self.queued_this_time = 0

        for _ in range (self.args.max_queue):
            smiles_list = self.read_smiles_chunk(standard_file,self.args.chunk_size)
            if len(smiles_list) > 0:
                self.process_queue.put(smiles_list)
                self.num_queued += 1
                self.queued_this_time +=1
                if len(smiles_list) < self.args.chunk_size:
                    smiles_to_process = False
                    break
            else:
                smiles_to_process = False
                break


        while (self.num_queued + self.num_requeued > self.num_processed):

            #Check if something has been processed
            while not self.results_queue.empty():

                node_list = self.results_queue.get()
                self.num_processed += 1
                # As replies are processed, room is created in the queue again
                if self.queued_this_time > 0:
                    self.queued_this_time -= 1

                # Process results
                for node_holder in node_list:
                    self.process_results(node_holder, max_frags=self.args.max_frag, verbosity=self.args.verbosity)

            #TODO Enough? - This acts as a limit but is not precise due to earlier chunking process
            if self.args.limit and self.smiles_read > self.args.limit:
                break

            # If there are not too many items queued then refill from queue
            if (self.queued_this_time < self.args.max_queue) and smiles_to_process:
                for _ in range(self.args.max_queue-self.queued_this_time):
                    smiles_list = self.read_smiles_chunk(standard_file,self.args.chunk_size)
                    if len(smiles_list) > 0:
                        self.process_queue.put(smiles_list)
                        self.num_queued += 1
                        if len(smiles_list) < self.args.chunk_size:
                            smiles_to_process = False
                            break
                    else:
                        smiles_to_process = False
                        break
                self.queued_this_time = 0

            #TODO re-enable this somehow - problem is that the code loops waiting for requests/replies
            #TODO so it can remain on one num_processed for a while.
            #if self.args.report_interval > 0 and self.num_processed % self.args.report_interval == 0:
            #    print("Number Queued {}, Queue Processed {}".format (self.num_queued, self.num_processed))

        #EndWhile

        print("Number Requests {}, Number Requeued {}, Number Responses {}"
              .format(self.num_queued, self.num_requeued, self.num_processed))
        print('Fragmentation Control Thread End')

        return self.num_processed

class FileWriter:
    """Class FileWriter
    Purpose:
        Contains all the output file writing for the fragmenting process.

    Methods:

    """

    node_count = 0
    edge_count = 0
    reject_count = 0

    def __init__(self, args, nodes_f, edges_f, rejects_f) -> None:
        '''Initialises the object.'''
        self.args = args
        self.nodes_f = nodes_f
        self.edges_f = edges_f
        self.rejects_f = rejects_f


    def get_node_count(self) -> int:
        return self.node_count

    def get_edge_count(self) -> int:
        return self.edge_count

    def get_reject_count(self) -> int:
        return self.reject_count

    def write_node(self, node, time_ms, num_children, num_edges):
        """Write to Note File.

        """
        # print("writing node", node.SMILES)
        self.node_count +=1
        self.nodes_f.write(','.join([node.SMILES, str(node.HAC), str(node.RAC), str(num_children), str(num_edges), str(time_ms)]) + '\n')

    def write_edge(self, edge):
        """Write to Edge File.

        """
        # print("writing edge", edge.get_label())
        self.edge_count +=1
        self.edges_f.write(edge.as_csv() + '\n')

    def write_reject(self, smiles):
        """Write to Reject File.

        """
        self.reject_count +=1
        self.rejects_f.write(smiles + '\n')

