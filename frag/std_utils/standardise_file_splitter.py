#!/usr/bin/env python

"""standardise_file_splitter.py

Splits a standardised input file into smaller files of the same type.
If a token is supplied then split is made at token over a number of lines, otherwise it is by line.

Example usages:

By token:
python -m frag.std_utils.standardise_file_splitter data/sdf/test.sdf.gz --file_prefix chunk_ chunk_size output_dir \
--token $$$$

By line:
python -m frag.std_utils.standardise_file_splitter data/xchem/dsip/dsip.txt.gz --file_prefix chunk_ chunk_size \
output_dir --header

Duncan Peacock December 2020 Initial Version

"""

import argparse
import gzip
import logging
import os
import sys

# Configure basic logging
logger = logging.getLogger('standardise_sdf')
out_hdlr = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s %(levelname)s # %(message)s',
                              '%Y-%m-%dT%H:%M:%S')
out_hdlr.setFormatter(formatter)
out_hdlr.setLevel(logging.INFO)
logger.addHandler(out_hdlr)
logger.setLevel(logging.INFO)


def error(msg):
    """Prints an error message and exists.

    :param msg: The message to print
    """
    logger.error('ERROR: {}'.format(msg))
    sys.exit(1)


def split_file(args):
    """Process the given file and split into chunks
    params:
        args from main.
    """

    num_processed = 0
    num_chunks = 1

    logger.info('Chunking %s...', args.input_file)
    logger.info(args)

    if args.input_file.endswith('.gz'):
        file = gzip.open(args.input_file,'rt')
    else:
        file = open(args.input_file,'r')

    start_chunk = 1
    output_filename = args.chunk_prefix + str(start_chunk).zfill(10)
    existing_file = os.path.join(args.output_dir, output_filename)

    while os.path.isfile(existing_file):
        start_chunk += 1
        output_filename = args.chunk_prefix + str(start_chunk).zfill(10)
        existing_file = os.path.join(args.output_dir, output_filename)

    output_file = open(os.path.join(args.output_dir, output_filename), 'wt')

    header = ''
    if args.header:
        header = file.readline()
        lines = file.readlines()[1:]
        output_file.write(header)
    else:
        lines = file.readlines()

    chunk_recs = 0
    new_file = False

    for line in lines:

        if new_file:
            start_chunk += 1
            output_filename = args.chunk_prefix + str(start_chunk).zfill(10)
            output_file = open(os.path.join(args.output_dir, output_filename), 'wt')
            if args.header:
                output_file.write(header)
            num_chunks += 1
            new_file = False

        output_file.write(line)

        if args.token == 'None':
            chunk_recs += 1
        else:
            # Check to see if the line is equal to a token before incrementing record counter
            if args.token in line:
                chunk_recs += 1

        if chunk_recs > args.chunk_size-1:
            new_file = True
            chunk_recs = 0

        num_processed += 1

    return num_processed, num_chunks


if __name__ == '__main__':

    parser = argparse.ArgumentParser('Standardise Input File Splitter')
    # python -m frag.std_utils.standardise_file_splitter data/sdf/test.sdf.gz --chunk_prefix chunk_ chunk_size /output
    # --token $$$$
    parser.add_argument('input_file',
                        help='The file to be split,')
    parser.add_argument('chunk_prefix',
                        help='The prefix to be used for the output files,'
                             ' will be appended by a number e.g. chunk_00001, chunk_00002 etc'
                             ' Note that its clever enough to continue if chunks exist ')
    parser.add_argument('chunk_size',
                        type=int,
                        help='The number of records to split at')
    parser.add_argument('output_dir',
                        help='The output directory')
    parser.add_argument('--token',
                        type=str, default='None',
                        help='If this is specified then the split is made at the given token rather than'
                             'by number of lines')
    parser.add_argument('--header',
                        action='store_true', default=False,
                        help='If this is set then each chunk will have the first line of the file written to it')
    args = parser.parse_args()

    if not os.path.isdir(args.output_dir):
        os.mkdir(args.output_dir)

    num_processed = 0
    num_chunks = 0
    validated = True

    if not os.path.isfile(args.input_file):
        print('No input file found')
        validated = False

    if args.chunk_size == 0:
        print('A chunksize must be specified')
        validated = False

    if validated:
        num_processed, num_chunks = split_file(args)

    # Summary
    logger.info('{:,} lines processed'.format(num_processed))
    logger.info('{:,} number chunk files created'.format(num_chunks))
