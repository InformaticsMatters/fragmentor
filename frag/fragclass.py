#!/usr/bin/env python
# Purpose
#
# Based on build_db_from_smiles.py, this module defines classes used in the
# fragmentor module.
#
# Duncan Peacock
# February 2020

from multiprocessing import Process

class FragProcess(Process):
    """Class FragProcess

    Purpose:
	    Manages the fragmention of a process queue element (consisting of a number of smiles).
		Returns a NodeHolder object to the results queue.

    Parameters:

    Methods:

    """

    def __init__(self, process_queue, results_queue) -> None:
        '''Initialises the object.'''
        super(FragProcess, self).__init__()
        self.process_queue = process_queue
        self.results_queue = results_queue

    def run(self) -> object:
        """ fragment process.

        Parameters:
            A list of SMILES to process.

        Returns:
            A list of NodeHolder objects.

        """
        print('fragment process')


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

    def __init__(self, args, process_queue, results_queue, f_writer) -> None:
        '''Initialises the object.'''
        super(FragController, self).__init__()
        self.args = args
        self.process_queue = process_queue
        self.results_queue = results_queue
        self.f_writer = f_writer

    def run(self) -> int:
        """Fragmentation Control Thread.

        Returns:
            the number of smiles processed.

        """
        print('fragmentation control thread')
        print(self.args.chunk_size)
        self.f_writer.write()

class FileWriter:
    """Class FileWriter
    Purpose:
        Contains all the output file writing for the fragmenting process.

    Methods:

    """

    def __init__(self, args, nodes_f, edges_f) -> None:
        '''Initialises the object.'''
        self.args = args
        self.nodes_f = nodes_f
        self.edges_f = edges_f

    def write(self):
        """Fragmentation Control Thread.

        Returns:
            the number of smiles processed.

        """
        print('Write')

