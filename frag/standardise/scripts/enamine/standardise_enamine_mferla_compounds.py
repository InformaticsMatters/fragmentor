#!/usr/bin/env python

"""standardise_enamine_mferla_compounds.py

Processes Compounds produced by mferla, and generates a 'standard'
tab-separated output.

We create a 'real-standardised-compounds.tab' file that contains a 1st-line
'header' formed from the _OUTPUT_COLUMNS list.

Alan Christie
September 2024
"""

import argparse
import glob
import logging
import os
import sys

from rdkit import RDLogger

from frag.utils import standardise_utils
from frag.std_utils import parser

# Configure basic logging
logger = logging.getLogger('real')
out_hdlr = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s %(levelname)s # %(message)s',
                              '%Y-%m-%dT%H:%M:%S')
out_hdlr.setFormatter(formatter)
out_hdlr.setLevel(logging.INFO)
logger.addHandler(out_hdlr)
logger.setLevel(logging.INFO)

# The columns in our output file.
# In this file we don't add any of our own.
_OUTPUT_COLUMNS = parser.STANDARD_COLUMNS

# The minimum number of columns in the input files and
# and a map of expected column names indexed by (0-based) column number.
#
# The 'standardised' files contain at least 2 columns...
#
# smiles    0
# idnumber  1

expected_min_num_cols = 2
smiles_col = 0
compound_col = 1
expected_input_cols = {compound_col: 'identifier',
                       smiles_col: 'smiles'}

# The output file.
output_filename = 'standardised-compounds.tab'

# The prefix we use in our fragment file
real_prefix = 'REAL:'

# The set of all original SMILES we've seen.
vendor_osmiles = set()
# A map of all the original SMILES indexed by their compound ID.
# Used to detect duplicate compounds.
vendor_compound_ids = {}

# Various diagnostic counts
num_vendor_mols = 0
num_vendor_molecule_failures = 0

# The line rate at which the process writes updates to stdout.
report_rate = 250000


def error(msg):
    """Prints an error message and exists.

    :param msg: The message to print
    """
    logger.error('ERROR: {}'.format(msg))
    sys.exit(1)


def standardise_vendor_compounds(output_file, file_name, limit):
    """Process the given file and standardise the vendor
    information, writing it as tab-separated fields to the output.

    As we load the Vendor compounds we 'standardise' the SMILES and
    determine whether they represent an isomer or not.

    :param output_file: The tab-separated standardised output file
    :param file_name: The (compressed) file to process
    :param limit: Limit processing to this number of values (or all if 0)
    :returns: The number of items processed
    """
    global vendor_compounds
    global num_vendor_mols
    global num_vendor_molecule_failures

    logger.info('Standardising %s...', file_name)

    line_num = 0
    num_processed = 0
    with open(file_name, 'rt', encoding='utf8') as input_file:

        # Check first line (a space-delimited header).
        # This is a basic sanity-check to make sure the important column
        # names are what we expect.

        hdr = input_file.readline()
        field_names = hdr.split()
        # Expected minimum number of columns...
        if len(field_names) < expected_min_num_cols:
            error('expected at least {} columns found {}'.
                  format(expected_input_cols, len(field_names)))
        # Check salient columns (ignoring case)...
        for col_num in expected_input_cols:
            actual_name = field_names[col_num].strip().lower()
            if actual_name != expected_input_cols[col_num]:
                error('expected "{}" in column {} found "{}"'.
                      format(expected_input_cols[col_num],
                             col_num,
                             actual_name))

        # Columns look right...

        for line in input_file:

            line_num += 1
            fields = line.split('\t')
            if len(fields) <= 1:
                continue

            if line_num % report_rate == 0:
                logger.info('...at compound {:,}'.format(line_num))

            osmiles = fields[smiles_col].split()[0]
            cid = fields[compound_col]
            compound_id = real_prefix + cid

            if osmiles in vendor_osmiles:
                # Potential duplicate?
                # A problem if the ID is already used and it's not the same SMILES
                if cid in vendor_compound_ids:
                    if vendor_compound_ids[cid] == osmiles:
                        logger.warning('Skipping duplicate compound %s %s', osmiles, cid)
                        continue
                    logger.warning('Found ID used for different compounds %s', cid)
                else:
                    logger.warning('Found duplicate compound with different IDs %s', osmiles)

                # If we get here the compound's been seen but it has the same ID
                # so we can skip it!

            vendor_osmiles.add(osmiles)
            vendor_compound_ids[cid] = osmiles

            # Standardise and update global maps...
            # And try and handle and report any catastrophic errors
            # from dependent modules/functions.

            std_info = standardise_utils.standardise(osmiles)
            if not std_info.std:
                num_vendor_molecule_failures += 1
                continue
            num_vendor_mols += 1

            # Write the standardised data

            output = [osmiles,
                      std_info.iso,
                      std_info.noniso,
                      std_info.hac,
                      compound_id]

            output_file.write('\t'.join(output) + '\n')

            # Enough?
            num_processed += 1
            if limit and num_processed >= limit:
                break

    return num_processed


if __name__ == '__main__':

    parser = argparse.ArgumentParser('Vendor Compound Standardiser (Enamine)')
    parser.add_argument('vendor_dir',
                        help='The Enamine vendor directory,'
                             ' containing the ".gz" files to be processed.')
    parser.add_argument('vendor_prefix',
                        help='The Enamine vendor file prefix,'
                             ' i.e. "June2018". Only files with this prefix'
                             ' in the vendor directory will be processed')
    parser.add_argument('output',
                        help='The output directory')
    parser.add_argument('--output-is-prefix',
                        action='store_true',
                        help='Use the output as filename prefix rather than'
                             ' a directory. This is useful in nextflow'
                             ' workflows')
    parser.add_argument('-l', '--limit',
                        type=int, default=0,
                        help='Limit processing to the first N molecules,'
                             ' process all otherwise')

    args = parser.parse_args()

    # Output is either s fixed name in an output directory
    # or a prefixed filename (without an output directory)
    if args.output_is_prefix:
        output_filename = f'{args.output}.{output_filename}'
    else:
        # Create the output directory
        if os.path.exists(args.output):
            logger.error('Output exists')
            sys.exit(1)
        os.mkdir(args.output)
        os.chmod(args.output, 0o777)
        output_filename = os.path.join(args.output, f'{output_filename}')

    # Suppress basic RDKit logging...
    RDLogger.logger().setLevel(RDLogger.ERROR)

    # Report any limiting...?
    if args.limit:
        logger.warning('Limiting processing to first {:,} molecules'.format(args.limit))

    # Before we open the output file
    # get a lit of all the input files (the prefix may be the same)
    # so we don't want our file in the list of files to be processed)
    real_files = glob.glob('{}/{}*'.format(args.vendor_dir,
                                              args.vendor_prefix))

    # Open the file we'll write the standardised data set to.
    # A text, tab-separated file.
    logger.info('Writing %s...', output_filename)
    num_processed = 0
    with open(output_filename, 'wt', encoding='utf8') as output_file:

        # Write the header...
        output_file.write('\t'.join(_OUTPUT_COLUMNS) + '\n')

        # Process all the Vendor files...
        for real_file in real_files:
            num_processed += standardise_vendor_compounds(output_file,
                                                          real_file,
                                                          args.limit)
            if args.limit and num_processed >= args.limit:
                break

    # Summary
    logger.info('{:,} vendor molecules'.format(num_vendor_mols))
    logger.info('{:,} vendor molecule failures'.format(num_vendor_molecule_failures))
