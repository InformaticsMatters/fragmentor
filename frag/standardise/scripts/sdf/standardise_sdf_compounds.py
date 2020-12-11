#!/usr/bin/env python

"""standardise_sdf_compounds.py

Processes a SDF file and generates a 'standard' tab-separated output.

We create a 'standardised-compounds.tab' file that contains a 1st-line
'header' formed from the _OUTPUT_COLUMNS list.

A field that contains the compound ID can be specified with the --id-field parameter.
If not specified the title line is used.
A prefix for the ID that is generated must be specified using the --prefix parameter.
e.g. is you specify --prefix FOO then the IDs are generated as FOO:<id-field-value>.

Example usage:
python -m frag.standardise.scripts.sdf.standardise_sdf_compounds --prefix FOO --id-field mr_id data/sdf Kinase_inhibs.sdf.gz foo

Tim Dudgeon November 2020 Initial Version

"""

import argparse
import glob
import gzip
import logging
import os
import sys
import traceback

from rdkit import Chem, RDLogger

from frag.std_utils import parser
from frag.utils import standardise_utils

# Configure basic logging
logger = logging.getLogger('standardise_sdf')
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

# The output file.
output_filename = 'standardised-compounds.tab'

# All the vendor compound IDs
vendor_compounds = set()
# A map of duplicate compounds and the number of duplicates.
# The index uses the vendor's original ID value, not our prefixed value.
duplicate_suffix = '-duplicate-'
vendor_duplicates = {}

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


def standardise_vendor_compounds(output_file, file_name, id_field, prefix, limit):
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
    global vendor_duplicates
    global duplicate_suffix
    global num_vendor_mols
    global num_vendor_molecule_failures

    logger.info('Standardising %s...', file_name)

    # If no field specified we use the title line which RDKit stores as the _Name field
    if not id_field:
        id_field = '_Name'

    line_num = 0
    num_processed = 0
    if file_name.endswith('.gz'):
        gz = gzip.open(file_name)
        supplr = Chem.ForwardSDMolSupplier(gz)
    else:
        supplr = Chem.ForwardSDMolSupplier(file_name)
    for mol in supplr:

        if not mol:
            # RDKit could not handle the record
            num_vendor_molecule_failures += 1
            continue
        try:
            osmiles = Chem.MolToSmiles(mol)
            vendor_id = mol.GetProp(id_field)
            compound_id = prefix + ':' + vendor_id

            # Add the compound (expected to be unique)
            # to our set of 'all compounds'.
            if compound_id in vendor_compounds:
                # Get the number of duplicates (default of 1)
                # using the vendor's original ID as a key
                duplicate_count = vendor_duplicates.setdefault(vendor_id, 1)
                compound_id += '{}{}'.format(duplicate_suffix, duplicate_count)
                # Increment for any further duplicates
                vendor_duplicates[vendor_id] = duplicate_count + 1
            else:
                vendor_compounds.add(compound_id)

            # Standardise and update global maps...
            # And try and handle and report any catastrophic errors
            # from dependent modules/functions.

            std_info = standardise_utils.standardise_mol(mol, osmiles)
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

        except:
            num_vendor_molecule_failures += 1
            traceback.print_exc()

        # Enough?
        num_processed += 1
        if limit and num_processed >= limit:
            break

    return num_processed


if __name__ == '__main__':

    parser = argparse.ArgumentParser('Vendor Compound Standardiser (SDF)')
    parser.add_argument('vendor_dir',
                        help='The SDF vendor directory,'
                             ' containing the ".gz" files to be processed.')
    parser.add_argument('vendor_prefix',
                        help='The SDF vendor file prefix,'
                             ' i.e. "iis_smiles". Only files with this prefix'
                             ' in the vendor directory will be processed')

    parser.add_argument('output',
                        help='The output directory')
    parser.add_argument('--id-field',
                        help='Name of the field for the compound ID. If not specified the title line is used')
    parser.add_argument('--prefix', required=True,
                        help='Prefix for the compound ID')
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
        output_filename = '{}.{}'.format(args.output, output_filename)
    else:
        # Create the output directory
        if os.path.exists(args.output):
            logger.error('Output exists')
            sys.exit(1)
        os.mkdir(args.output)
        os.chmod(args.output, 0o777)
        output_filename = os.path.join(args.output,
                                       '{}'.format(output_filename))

    # Suppress basic RDKit logging...
    RDLogger.logger().setLevel(RDLogger.ERROR)

    # Report any limiting...?
    if args.limit:
        logger.warning('Limiting processing to first {:,} molecules'.format(args.limit))

    sdf_files = glob.glob('{}/{}*'.format(args.vendor_dir,
                                                 args.vendor_prefix))

    # Open the file we'll write the standardised data set to.
    # A text, tab-separated file.
    logger.info('Writing %s...', output_filename)

    num_processed = 0
    with open(output_filename, 'wt') as output_file:

        # Write the header...
        output_file.write('\t'.join(_OUTPUT_COLUMNS) + '\n')

        # Process all the Vendor files...
        for sdf_file in sdf_files:
            num_processed += standardise_vendor_compounds(output_file,
                                                          sdf_file,
                                                          args.id_field,
                                                          args.prefix,
                                                          args.limit)
            if args.limit and num_processed >= args.limit:
                break

    # Summary
    logger.info('{:,} vendor molecules'.format(num_vendor_mols))
    logger.info('{:,} vendor molecule failures'.format(num_vendor_molecule_failures))
    logger.info('{:,} duplicate compounds'.format(len(vendor_duplicates)))
    for vendor_duplicate in vendor_duplicates:
        logger.info('Duplicate compound: {} x{}'.format(vendor_duplicate, vendor_duplicates[vendor_duplicate]))
