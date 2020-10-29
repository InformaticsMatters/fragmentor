#!/usr/bin/env python
# coding=utf-8

"""standardise_utils.py

Utils shared with the various standardise modules.

The module provides a method for rendering a vendor (original) SMILES
representation into a standard form, including attributes like the
heavy atom count. As each caller may augment the standard with their own
columns then they are required to construct the output file, and its
header.

The standard (mandatory) columns are available to all in the
STANDARD_COLUMNS list. This must form the first line of every file written
and they must be the first columns in the file.

Alan Christie
February 2019
"""

import logging
import hashlib
import os
from collections import namedtuple
from os import path

from rdkit import Chem
from frag.utils.rdkit_utils import standardize

# The tuple returned by calls to 'standardise()'.
# If the first field (std) is None then the standard cannot be used.
# All items are rendered as strings.
StandardInfo = namedtuple('StandardInfo', 'std iso noniso hac rac inchis inchik noniso_inchis noniso_inchik iso_inchis iso_inchik')

# Out logger
logger = logging.getLogger(__name__)

def standardise(osmiles):
    """Given a vendor (original) SMILES this method standardises
   it into a canonical form and returns a namedtuple that contains
   the standard form, the isomeric form, the non-isomeric form
   the heavy atom count (hac).

   :param osmiles: The original (non-standard) SMILES

   :return: A namedtuple containing the standard molecule
            representations and info. Errors are logged and, on error,
            the standard form will be returned as None.
   """
    mol = None
    try:
        mol = Chem.MolFromSmiles(osmiles)
    except Exception as e:
        logger.warning('MolFromSmiles(%s) exception: "%s"',
                       osmiles, e.message)
    return standardize(mol, osmiles)

def standardise(mol, osmiles):
    """Standaridise a RDKit mol object. A SMILES string that represents the molecule can also be passed in to allow
    better logging of errors.
    """
    global logger

    # Standardise and update global maps...
    # And try and handle and report any catastrophic errors
    # from dependent modules/functions.

    std = None
    iso = None
    noniso = None
    inchis = None
    inchik = None
    noniso_inchis = None
    noniso_inchik = None
    iso_inchis = None
    iso_inchik = None

    hac = 0
    rac = 0

    if not mol:
        logger.error('Got nothing from MolFromSmiles(%s).'
                     ' Skipping this Vendor compound', osmiles)

    if mol:

        # Got a molecule.
        #
        # Get the HAC and try to (safely) standardise,
        # and create isomeric an non-isomeric representations.

        hac = mol.GetNumHeavyAtoms()
        rac = 0
        for atom in mol.GetAtoms():
            if atom.IsInRing():
                rac += 1

        try:
            std = standardize(mol)
        except Exception as e:
            logger.warning('standardize(%s) exception: "%s"',
                           osmiles, e.message)
        if not std:
            logger.error('Got nothing from standardize(%s).'
                         ' Skipping this Vendor compound', osmiles)

    if std:

        # We have a standard representation,
        # Try to generate the isomeric version...

        try:
            iso = Chem.MolToSmiles(std, isomericSmiles=True, canonical=True)
        except Exception as e:
            logger.warning('MolToSmiles(%s, iso) exception: "%s"',
                           osmiles, e.message)
        if not iso:
            logger.error('Got nothing from MolToSmiles(%s, iso).'
                         ' Skipping this Vendor compound', osmiles)


    if std:

        # We have a standard representation,
        # Try to generate the non-isomeric version...

        try:
            noniso = Chem.MolToSmiles(std, isomericSmiles=False, canonical=True)
            nonisoMol = Chem.MolFromSmiles(noniso)
            try:
                inchis, inchik = gen_inchi(nonisoMol, '')
                iso_inchis, iso_inchik = gen_inchi(std, '/FixedH')
                noniso_inchis, noniso_inchik = gen_inchi(nonisoMol, '/FixedH')
            except Exception as e:
                logger.warning('gen_inchi exception for %s', noniso)

        except Exception as e:
            logger.warning('MolToSmiles(%s, noniso) exception: "%s"',
                           osmiles, e.message)
        if not noniso:
            logger.error('Got nothing from MolToSmiles(%s, noniso).'
                         ' Skipping this Vendor compound', osmiles)


    # If anything went wrong, set std to None.
    # It's "all-or-nothing"...
    if not iso or not noniso:
        std = None

    return StandardInfo(std, iso, noniso, str(hac), str(rac), inchis, inchik, noniso_inchis, noniso_inchik, iso_inchis, iso_inchik)

def gen_inchi(mol, opts):
    inchis = Chem.inchi.MolToInchi(mol, opts)
    inchik = Chem.inchi.InchiToInchiKey(inchis)
    return inchis, inchik

def write_data(tree_root, vendor, osmiles, std_info, compound_id, vendor_paths_file):

    global logger

    new_inchis = 0
    new_noniso = 0
    new_iso = 0

    inchis = std_info.inchis
    inchik = std_info.inchik
    prefix = inchik[:2]
    level1 = "/".join([tree_root, prefix]) # first 2 chars of InChiKey
    level2 = "/".join([level1, inchik]) # full InChiKey
    md5_inchis = hashlib.md5(inchis.encode('utf-8'))
    md5_noniso = hashlib.md5(std_info.noniso.encode('utf-8'))
    level3 = "/".join([level2, md5_inchis.hexdigest()]) # hash of InChi string
    inchidata = "/".join([level3, 'inchi'])
    level4 = "/".join([level3, md5_noniso.hexdigest()]) # hash of non-isomeric smiles
    nonisodata = "/".join([level4, 'moldata.yml'])



    if not path.exists(level1): # first 2 chars of InChiKey
        #logger.info('Creating dir at level 1 %s', level1)
        os.mkdir(level1)
    # else:
    #     logger.info('Level 1 already exists %s', level1)

    if not path.exists(level2): # full InChiKey
        #logger.info('Creating dir at level 2 %s', level2)
        os.mkdir(level2)
    # else:
    #     logger.info('Level 2 already exists %s', level2)

    if not path.exists(level3): # hash of InChi string
        #logger.info('Creating dir at level 3 %s', level3)
        os.mkdir(level3)
        f = open(inchidata, "wt")
        f.write(inchis)
        f.close()
        new_inchis += 1
    # else:
    #     logger.info('Level 3 already exists %s', level3)

    if not path.exists(level4): # hash of non-isomeric smiles
        #logger.info('Creating dir at level 4 %s', level4)
        os.mkdir(level4)
        f = open(nonisodata, "wt")
        f.write("smiles: {}\n".format(std_info.noniso))
        f.write("inchi: {}\n".format(std_info.noniso_inchis))
        f.write("inchikey: {}\n".format(std_info.noniso_inchik))
        f.write("hac: {}\n".format(std_info.hac))
        f.close()
        new_noniso += 1
    # else:
    #     logger.info('Level 4 already exists %s', level4)

    if (std_info.noniso == std_info.iso):
        vendorlevel = level4
    else:
        md5_iso = hashlib.md5(std_info.iso.encode('utf-8'))
        level5 = "/".join([level4, md5_iso.hexdigest()])
        vendorlevel = level5
        if not path.exists(level5): # hash of isomeric smiles
            isodata = "/".join([level5, 'moldata.yml'])
            os.mkdir(level5)
            f = open(isodata, "wt")
            f.write("smiles: {}\n".format(std_info.iso))
            f.write("inchi: {}\n".format(std_info.iso_inchis))
            f.write("inchikey: {}\n".format(std_info.iso_inchik))
            f.close()
            new_iso += 1
        # else:
        #     logger.info('Level 5 already exists %s', level5)

    vendordata = "/".join([vendorlevel, vendor + '.yml'])

    # if path.exists(vendordata): # vendor data
    #     logger.info('Appending to existing data for %s', vendordata)
    # else:
    trimmed = vendorlevel[len(tree_root)+1:]
    vendor_paths_file.write(trimmed + "\n")

    f = open(vendordata, "at")
    f.write("id: {}\nsmiles: {}\n".format(compound_id, osmiles))
    f.close()

    return new_inchis, new_noniso, new_iso
