#!/usr/bin/env python
# Purpose
#
# Based on build_db_from_smiles.py, this module defines classes used in the
# fragmentor module.
#
# Duncan Peacock
# February 2020

import os
import sys
import collections
from multiprocessing import Process
from frag.network.models import NodeHolder, Attr
from frag.utils.network_utils import build_network

class FragData():
    """Class FragData

    Purpose:

    """

    def __init__(self):
        self.parent_data = collections.OrderedDict()
        self.nodes_map = {}

    def add_nodes(self, nodes, input_smiles):

        for node in nodes:
            self.nodes_map[node.SMILES] = node
        for smiles in input_smiles:
            node = self.nodes_map[smiles]
            self.parent_data[smiles] = ParentData(smiles, node.HAC, node.RAC)

    def add_edges(self, edge):

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
    """Class ParentData

    Purpose:

    """

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
           Fragdata object with Node/Edge data to write to files.

        """
        # Note that in this version, only one SMILES is sent in here.
        # There seemed to be some strange issues with edges if a combined node holder is used TBI
        attrs = []
        attr = Attr(smiles, ["EM"])
        attrs.append(attr)

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
                frag_data = FragData()
                for smiles in chunk_of_smiles:
                    node_holder = self.fragment_mol(smiles, verbosity=self.args.verbosity)
                    frag_data.add_nodes(node_holder.get_nodes(), [smiles])
                    for edge in node_holder.get_edges():
                        frag_data.add_edges(edge)
                self.results_queue.put(frag_data)


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
    smiles_requeued = 0
    max_smiles_requeued = 0
    # This is just an indication for sizing - divide by num_processed.
    total_size_of_result = 0
    max_size_of_result = 0
    cache = set()
    reprocess_buffer = set()
    already_processed = 0

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

    def write_data(self, parent_data):
        """Write_data

        node and edge data
        return: the set of smiles that needs further processing.

        """

        children_to_process = set()
        p_smiles = parent_data.smiles
        num_children = len(parent_data.children)
        num_edges = parent_data.edge_count

        # if no edges then this is a leaf node with no children so we just write the node
        if num_edges == 0:
            if p_smiles not in self.cache:
                self.f_writer.write_node(p_smiles, parent_data.hac, parent_data.rac, 0, num_edges)
                self.cache.add(p_smiles)
            return children_to_process
        else:
            # so we have edges to process
            if p_smiles not in self.cache:
                self.f_writer.write_node(p_smiles, parent_data.hac, parent_data.rac, 0, num_edges)
                self.cache.add(p_smiles)

                for c_smiles in parent_data.children:
                    labels = parent_data.children[c_smiles]
                    for label in labels:
                        self.f_writer.write_edge(p_smiles, c_smiles, label)
                    if c_smiles not in self.cache:
                        children_to_process.add(c_smiles)

            return children_to_process

    def write_process_queue (self, buffer, main_queue):
        """Write a chunk to the process queue and set flags.

        """
        self.process_queue.put(buffer)
        self.queued_this_time += 1
        if main_queue:
            self.num_queued += 1
        else:
            self.num_requeued += 1

        if self.args.report_interval > 0 and \
            (self.num_queued + self.num_requeued) % self.args.report_interval == 0:
            print("Number Queued {}, Number Requeued {}".format(self.num_queued, self.num_requeued))

    def reprocess_smiles_chunk(self, smiles_to_process, chunk_size, flush):
        """Create a chunk of SMILES to reprocess by buffering.
           If there are no SMILES in queue, then flush the requeue buffer.

        """
        buffer_size = len(self.reprocess_buffer)
        if flush and buffer_size > 0:
            self.write_process_queue (self.reprocess_buffer, False)
            self.reprocess_buffer.clear()
        else:
            reprocess_count = len(smiles_to_process)
            if reprocess_count == 0:
                return

            self.smiles_requeued += reprocess_count
            #statistics
            if reprocess_count > self.max_smiles_requeued:
                self.max_smiles_requeued = reprocess_count

            # small performance improvement
            if (buffer_size + reprocess_count) < chunk_size:
                self.reprocess_buffer = self.reprocess_buffer.union(smiles_to_process)
                buffer_size += reprocess_count
                return

            # write buffer on changeover
            for smiles in smiles_to_process:
                if buffer_size < chunk_size:
                    self.reprocess_buffer.add(smiles)
                    buffer_size += 1
                else:
                    self.write_process_queue(self.reprocess_buffer, False)
                    self.reprocess_buffer.clear()
                    self.reprocess_buffer.add(smiles)
                    buffer_size = 1

    def process_results(self, frag_data, max_frags=0, verbosity=0):
        """Processes a NodeHolder object from the results queue.

        """

        # Set of ALL parent smiles including ones that are in the cache
        data = frag_data.get_parent_data_list()

        # Loop through remaining parent_data not in cache.
        for parent_data in data:
            if parent_data.smiles in self.cache:
                continue

            num_children = len(parent_data.children)
            if 0 < max_frags < num_children:
               self.f_writer.write_reject(parent_data.smiles)
               continue

            # This writes the parent nodes to the cache.
            children_needing_processing = self.write_data(parent_data)
            self.reprocess_smiles_chunk(children_needing_processing, self.args.chunk_size, False)

    def read_smiles_chunk(self, standard_file, no_of_lines, process) -> bool:
        """Read a chunk of SMILES from the standard file (and write to the process queue.

        Returns:
            A list of SMILES to add to the process queue.

        """
        smiles_list = []
        more_smiles = True

        while True and len(smiles_list) < no_of_lines:
            # read a single line
            line = standard_file.readline().strip()
            if line:
                if line not in self.cache:
                    smiles_list.append(line)
            else:
                more_smiles = False
                break

        no_of_smiles = len(smiles_list)
        if no_of_smiles == 0:
            return False

        self.smiles_read += no_of_smiles
        if process:
            self.write_process_queue(smiles_list, True)
        return more_smiles


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
            smiles_to_process = self.read_smiles_chunk(standard_file, self.args.skip, False).size()
        self.queued_this_time = 0

        for _ in range (self.args.max_queue):
            smiles_to_process = self.read_smiles_chunk(standard_file,self.args.chunk_size, True)
            # If reached limit then stop.
            if self.args.limit and self.smiles_read > self.args.limit:
                smiles_to_process = False
                break

        while (self.num_queued + self.num_requeued > self.num_processed):

            #Check if something has been processed
            while not self.results_queue.empty():

                #node_list = self.results_queue.get()
                frag_data = self.results_queue.get()
                self.num_processed += 1
                # As replies are processed, room is created in the queue again
                if self.queued_this_time > 0:
                    self.queued_this_time -= 1
                size_of_result =  sys.getsizeof(frag_data)
                self.total_size_of_result += size_of_result
                if size_of_result > self.max_size_of_result:
                    self.max_size_of_result = size_of_result

                self.process_results(frag_data, max_frags=self.args.max_frag, verbosity=self.args.verbosity)

            # Flushes remaining smiles in reprocess buffer
            self.reprocess_smiles_chunk([], self.args.chunk_size, True)

            if self.args.limit and self.smiles_read > self.args.limit:
                smiles_to_process = False

            # If there are not too many items queued then refill from queue
            if (self.queued_this_time < self.args.max_queue) and smiles_to_process:
                for _ in range(self.args.max_queue-self.queued_this_time):
                    smiles_to_process = self.read_smiles_chunk(standard_file, self.args.chunk_size, True)
                self.queued_this_time = 0

        #EndWhile

        print("Number Requests {}, Number Requeued {}, Number Responses {}"
              .format(self.num_queued, self.num_requeued, self.num_processed))
        print("SMILES {}, SMILES Requeued {}, Max SMILES to requeued {}"
              .format(self.smiles_read, self.smiles_requeued, self.max_smiles_requeued))
        print("Average Results Size {}, Max Result Size {}"
              .format((self.total_size_of_result / self.num_processed), self.max_size_of_result))
        print("Found in Cache and removed from Reprocessing {}"
              .format(self.already_processed))

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

    def write_node(self, p_smiles, hac, rac, num_children, num_edges):
        """Write to Node File.
        """
        self.node_count +=1
        self.nodes_f.write(','.join([p_smiles, str(hac), str(rac), str(num_children), str(num_edges)]) + '\n')

    def write_edge(self, p_smiles, c_smiles, label):
        """Write to Edge File.
        """
        # print("writing edge", edge.get_label())
        self.edges_f.write(','.join([p_smiles, c_smiles, label]) + '\n')
        self.edge_count +=1

    def write_reject(self, smiles):
        """Write to Reject File.

        """
        self.reject_count +=1
        self.rejects_f.write(smiles + '\n')

