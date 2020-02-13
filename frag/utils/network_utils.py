from rdkit import Chem

try:
    from rdkit.Chem import fdMCS as MCS
except ImportError:
    from rdkit.Chem import MCS
from rdkit.Chem import AllChem, Draw
from tqdm import tqdm
import timeit
import os

SMARTS_PATTERN = "[*;R]-;!@[*]"

""" contribution from Hans de Winter """
from rdkit import Chem
from rdkit.Chem import AllChem

# Instrumentation to expose the 'retrieve', 'create' and 'total' durations,
# and the number of nodes and edges for each molecule processed
# by build_network().
#
# Data is flushed and synchronised to disc after each molecule
# (so that data is kept even if the process dies). So anticipate a
# drop in throughput if you enable this feature.
#
# Alan Christie, Sep 2018
ENABLE_BUILD_NETWORK_LOG = os.environ.get('ENABLE_BUILD_NETWORK_LOG')
BUILD_NETWORK_LOG = 'build-network.log'
CHILD_DEPTH = 0


def _InitialiseNeutralisationReactions():
    patts = (
        # Imidazoles
        ("[n+;H]", "n"),
        # Amines
        ("[N+;!H0]", "N"),
        # Carboxylic acids and alcohols
        ("[$([O-]);!$([O-][#7])]", "O"),
        # Thiols
        ("[S-;X1]", "S"),
        # Sulfonamides
        ("[$([N-;X2]S(=O)=O)]", "N"),
        # Enamines
        ("[$([N-;X2][C,N]=C)]", "N"),
        # Tetrazoles
        ("[n-]", "[nH]"),
        # Sulfoxides
        ("[$([S-]=O)]", "S"),
        # Amides
        ("[$([N-]C=O)]", "N"),
    )
    return [(Chem.MolFromSmarts(x), Chem.MolFromSmiles(y, False)) for x, y in patts]


_reactions = None


def NeutraliseCharges(mol, reactions=None, as_mol=False):
    global _reactions
    smiles = Chem.MolToSmiles(mol, True)
    if reactions is None:
        if _reactions is None:
            _reactions = _InitialiseNeutralisationReactions()
        reactions = _reactions
    replaced = False
    for i, (reactant, product) in enumerate(reactions):
        while mol.HasSubstructMatch(reactant):
            replaced = True
            rms = AllChem.ReplaceSubstructs(mol, reactant, product)
            mol = rms[0]
    if replaced:
        if as_mol:
            return (mol, True)
        else:
            return (Chem.MolToSmiles(mol, True), True)
    else:
        if as_mol:
            return (mol, True)
        else:
            return (smiles, False)


def conv_at_xe(x):
    return NeutraliseCharges(Chem.MolFromSmiles(x.replace("[At]", "[Xe]")))[0]


def write_results(input_dict):
    """
    Write out a structure list of results.
    :param input_dict:
    :return:
    """
    out_imgs = {}
    for mol in input_dict:
        mols = []
        for new_mols in input_dict[mol]:
            m = Chem.MolFromSmiles(new_mols)
            mols.append(m)
        if len(mols) > 2:
            p = Chem.MolFromSmarts(MCS.FindMCS(mols).smarts)
            AllChem.Compute2DCoords(p)
            for m in mols:
                AllChem.GenerateDepictionMatching2DStructure(m, p)
        out_imgs[mol] = Draw.MolsToGridImage(mols, useSVG=True)
        # Write out the image
    return out_imgs


def get_frag_list(str_find, input_mol):
    """
    Get the list of fragments
    :param str_find:
    :param input_mol:
    :return:
    """
    return [
        x.replace(str_find, "Xe")
        for x in Chem.MolToSmiles(input_mol, isomericSmiles=True).split(".")
    ]


def get_comb_index(bi_1, bi_2):
    """

    :param bi_1:
    :param bi_2:
    :return:
    """
    return bi_1 + (100 * (bi_2 + 1))


def ret_comb_index(bi_tot, get_indices=False, isotope=None):
    """

    :param bi_tot:
    :return:
    """
    bi_1 = int(str(bi_tot)[-2:])
    bi_2 = int(str(bi_tot)[0:-2])
    if get_indices:
        return (bi_1, bi_2 - 1, isotope)
    else:
        return (bi_1, bi_2 - 1)


def get_fragments(input_mol, iso_labels=True, get_index_iso_map=False):
    """
    Find the frgments for a given molecule
    :param input_mol:
    :return:
    """
    index_isotope_map = {}
    atom_indices = input_mol.GetSubstructMatches(Chem.MolFromSmarts(SMARTS_PATTERN))
    if atom_indices and iso_labels:
        counter = 100
        labels = []
        bs = []
        for bi in atom_indices:
            b = input_mol.GetBondBetweenAtoms(bi[0], bi[1])
            if counter in index_isotope_map:
                index_isotope_map[counter].append(b.GetIdx())
            else:
                index_isotope_map[counter] = [b.GetIdx()]
            labels.append((counter, counter))
            bs.append(b.GetIdx())
            counter += 1
        input_mol = Chem.FragmentOnBonds(input_mol, bs, dummyLabels=labels)
    elif atom_indices:
        bs = []
        labels = []
        for bi in atom_indices:
            b = input_mol.GetBondBetweenAtoms(bi[0], bi[1])
            bs.append(b.GetIdx())
            comb_index = get_comb_index(bi[0], bi[1])
            labels.append((comb_index, comb_index))
        input_mol = Chem.FragmentOnBonds(input_mol, bs, dummyLabels=labels)
        return get_frag_list(str_find="*", input_mol=input_mol)
    if get_index_iso_map:
        return get_frag_list(str_find="*", input_mol=input_mol), index_isotope_map
    else:
        return get_frag_list(str_find="*", input_mol=input_mol)


def get_num_ring_atoms(input_mol):
    """
    Get the number of ring atoms
    :param input_mol:
    :return:
    """
    num_ring_atoms = 0
    split_index = 0
    split_indices = []
    for atom in input_mol.GetAtoms():
        if atom.IsInRing():
            num_ring_atoms += 1
        else:
            split_indices.append(split_index)
        split_index += 1
    return num_ring_atoms, split_indices


def simplified_graph(input_smiles, iso_flag=True):
    """
    The simplified graph representation for the edge uses the daylight
function dt_molgraph to
1) set all bond orders to one, 2) remove aromaticity, 3) set the hydrogen count,
4) remove charges and set masses to zero. 5)  Additionally we set every ring atom element to carbon.
This is Anthony Bradley's version of the algorithm.
    :param input_smiles: the input smiles to simplify
    :param iso_flag: a boolean to determine if we should include chirality etc
    :return:
    """
    mol = Chem.MolFromSmiles(input_smiles)
    for atom in mol.GetAtoms():
        atom.SetFormalCharge(0)
        if atom.IsInRing():
            atom.SetAtomicNum(6)
            atom.SetIsAromatic(False)
        for bond in atom.GetBonds():
            bond.SetBondType(Chem.BondType.SINGLE)
        if not iso_flag:
            atom.SetChiralTag(Chem.ChiralType.CHI_OTHER)
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def recombine_edges(output_edges):
    """
    Recombine a list of edges based on their rules.
    Recombines identical Xe isotopes. Remove isotopes.
    :param output_edges:
    :return: A SMILES string or None on failure
    """
    mol = Chem.MolFromSmiles(".".join(output_edges))
    if mol is None:
        return None
    # Dictionary of atom's to bond together and delete if they come in pairs
    iso_dict = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 54:
            # Get the isotope
            iso = atom.GetIsotope()
            if iso in iso_dict:
                iso_dict[iso].append(get_info(atom))
            else:
                iso_dict[iso] = [get_info(atom)]
    mw = Chem.RWMol(mol)
    # Add bonds first
    del_indices = []
    for isotope in iso_dict:
        if len(iso_dict[isotope]) > 1:
            mw.AddBond(
                iso_dict[isotope][0][1], iso_dict[isotope][1][1], Chem.BondType.SINGLE
            )
            del_indices.append(iso_dict[isotope][0][0])
            del_indices.append(iso_dict[isotope][1][0])
    # Now delete atoms
    del_count = 0
    for atom_index in sorted(del_indices):
        mw.RemoveAtom(atom_index - del_count)
        del_count += 1
    Chem.SanitizeMol(mw)
    return Chem.MolToSmiles(mw, isomericSmiles=True)


def rebuild_smi(input_list, ring_ring):
    """
    Rebuild a SMILES
    :param input_list: the list of fragments to be rebuilt
    :param ring_ring: a boolean - indicating if this is a ring ring split
    :return: A SMILES String or None on failure
    """
    if ring_ring:
        rebuilt_smi = ".".join(input_list)
    else:
        rebuilt_smi = recombine_edges(input_list)
    return rebuilt_smi


def make_child_mol(rebuilt_smi):
    """
    Make the child molecule
    :param rebuilt_smi:
    :return:
    """
    mol = Chem.MolFromSmiles(rebuilt_smi)
    # TODO proper warning messages for thsee two exceptions
    if mol is None:
        return None
    # TODO - check that non isomeric is ok here
    new_smi = Chem.MolToSmiles(mol, isomericSmiles=False).replace("[Xe]", "[H]")
    mol = Chem.MolFromSmiles(new_smi)
    # TODO proper warning messages for thsee two exceptions
    if mol is None:
        return None
    return Chem.CanonSmiles(new_smi)


def get_info(atom):
    """
    Get the needed info for an atom
    :param atom:
    :return:
    """
    return [atom.GetIdx(), atom.GetNeighbors()[0].GetIdx()]


def get_type(smiles):
    """
    Get the Type of a given entity
    :param smiles: the input SMILES
    :return: the type (FG or RING)
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[*;R]")):
        return "RING"
    return "FG"


def get_driver(url="neo4j"):
    """
    Get the driver to the network connection
    :return: the driver for the graphdabase
    """
    # No auth on the database
    from neo4j.v1 import GraphDatabase

    driver = GraphDatabase.driver("bolt://" + url + ":7687")
    return driver


def get_ring_ring_splits(input_mol, labels=False, do_comb_index=False):
    """
    Get and break Atom-Atom pairs in two different rings.
    :param input_mol:
    :return:
    """
    # TODO Fix for fused e.g.s
    RI = input_mol.GetRingInfo()
    rings = RI.AtomRings()
    out_mols = []
    bonds = [item for sublist in RI.BondRings() for item in sublist]
    bs = []
    for bond in input_mol.GetBonds():
        if bond.GetIdx() in bonds:
            continue
        id_one = bond.GetBeginAtomIdx()
        id_two = bond.GetEndAtomIdx()
        # Now find all pairs that are in both
        for ring in rings:
            if id_one in ring:
                for ring_two in rings:
                    if ring == ring_two:
                        continue
                    if id_two in ring_two:
                        bs.append(bond)
    if bs:
        for b in bs:
            if labels:
                nm = Chem.FragmentOnBonds(
                    input_mol,
                    [b.GetIdx()],
                    dummyLabels=[(b.GetBeginAtomIdx(), b.GetEndAtomIdx())],
                )
                mols = [
                    x.replace("*", "Xe")
                    for x in Chem.MolToSmiles(nm, isomericSmiles=True).split(".")
                ]
            elif do_comb_index:
                print(b)

                comb_index = get_comb_index(b.GetBeginAtomIdx(), b.GetEndAtomIdx())
                nm = Chem.FragmentOnBonds(
                    input_mol, [b.GetIdx()], dummyLabels=[(comb_index, comb_index)]
                )
                mols = [
                    x.replace("*", "Xe")
                    for x in Chem.MolToSmiles(nm, isomericSmiles=True).split(".")
                ]
            else:
                nm = Chem.FragmentOnBonds(input_mol, [b.GetIdx()], dummyLabels=[(1, 1)])
                # Only takes first
                mols = [
                    x.replace("1*", "Xe")
                    for x in Chem.MolToSmiles(nm, isomericSmiles=True).split(".")
                ]
            out_mols.append(mols)

    return out_mols


def add_child_and_edge(new_list, input_node, excluded_smi, node_holder,
                       ring_ring=False, recurse=True):
    """
    :param input_pair:
    :return:
    """
    # Rebuild the molecule
    rebuilt_smi = rebuild_smi(new_list, ring_ring)
    if rebuilt_smi is None:
        return
    # Turn into child molecule
    child_smi = make_child_mol(rebuilt_smi)
    if child_smi is None:
        return
    # Now generate the edges with input and this node
    # (catering for node failure)
    node, is_new = node_holder.create_or_retrieve_node(child_smi)
    if node:
        node_holder.create_or_retrieve_edge(excluded_smi,
                                            rebuilt_smi,
                                            input_node,
                                            node)
        # Now generate the children for this too
        if is_new and recurse:
            create_children(node, node_holder)


def canon_input(smi, isomericSmiles=True):
    # Decharge in this step too
    iso_smiles = NeutraliseCharges(Chem.MolFromSmiles(smi))[0]
    if not isomericSmiles:
        return Chem.MolToSmiles(Chem.MolFromSmiles(iso_smiles), isomericSmiles=False)
    else:
        return iso_smiles


def create_children(input_node, node_holder, max_frag=0, smiles=None, log_file=None, recurse=True):
    """
    Create a series of edges from an input molecule. Iteratively
    :param input_node:
    :param max_frag: Max initial fragments (or no limit if 0)
    :param smiles: A SMILES string, for log/diagnostics only.
                   Written to the log if specified
    :param log_file: A file if information is to be logged, otherwise None
    :return: A tuple, the number of direct children and the number of ring-ring splits
    """
    fragments = get_fragments(input_node.RDMol)

    # Ditch processing if too few fragments
    # or (if a maximum is defined) too many.
    #
    # All exclusions are written to a log file (if provided)
    # using 'X' as a prefix followed by letters that indicate the
    # excluded reason (MF=MaxFrag) and then a period
    # followed by a number (in this case the fragment count)
    num_fragments = len(fragments)
    if num_fragments < 2:
        return num_fragments, 0
    elif max_frag > 0 and num_fragments > max_frag:
        if log_file:
            log_file.write('XMF.%s,%s\n' % (num_fragments, smiles))
        return num_fragments, 0

    # OK if we get here.
    # If we have a log file write the fragment count.
    if log_file:
        log_file.write('F.%s,%s\n' % (num_fragments, smiles))

    # Get all ring-ring splits
    ring_ring_splits = get_ring_ring_splits(input_node.RDMol)

    num_ring_ring_splits = len(ring_ring_splits)
    num_total = num_ring_ring_splits + num_fragments
    if num_total < 2:
        return num_fragments, num_ring_ring_splits

    # Now remove one item on each iteration
    for i in range(num_fragments):
        new_list = []
        for j, item in enumerate(fragments):
            if i == j:
                excluded_smi = item
                continue
            new_list.append(item)
        add_child_and_edge(new_list, input_node, excluded_smi, node_holder, recurse=recurse)

    if ring_ring_splits:
        for ring_ring_split in ring_ring_splits:
            add_child_and_edge(
                ring_ring_split, input_node, "[Xe]", node_holder, ring_ring=True, recurse=recurse
            )

    return num_fragments, num_ring_ring_splits

def neutralise_3d_mol(input_mol):
    neutral_mol = NeutraliseCharges(Chem.MolFromMolBlock(input_mol), as_mol=True)[0]
    return AllChem.ConstrainedEmbed(neutral_mol, Chem.MolFromMolBlock(input_mol))


def write_data(output_dir, node_holder, attrs):
    """Write the node holder (nodes and edges) and, if attrs is not None
    the attributes.

    :param output_dir: The output directory
    :param node_holder: Written as nodes.txt and edges.txt
    :param attrs: Written as attributes.txt if not None
    """
    out_f = open(os.path.join(output_dir, "nodes.txt"), "w")
    for node in node_holder.node_list:
        out_f.write(str(node))
        out_f.write("\n")
    out_f.close()
    out_f = open(os.path.join(output_dir, "edges.txt"), "w")
    for edge in node_holder.get_edges():
        out_f.write(str(edge))
        out_f.write("\n")
    out_f.close()
    if attrs:
        out_f = open(os.path.join(output_dir, "attributes.txt"), "w")
        for attr in attrs:
            out_f.write(str(attr))
            out_f.write("\n")
        out_f.close()


def write_data_as_csv(output_dir, node_holder):
    """Write the node holder (nodes and edges). The output files are
    header-less CSV files. A final column for the node namespace and a
    label for the edge is added so that written nodes have 6 columns
    (2 are added for an initially empty compound ID and a label of 'F2')
    and edges have 4 (one added for a label of 'FRAG').

    :param output_dir: The output directory
    :param node_holder: Written as nodes.txt and edges.txt
    """
    out_f = open(os.path.join(output_dir, "nodes.csv"), "w")
    for node in node_holder.node_list:
        out_f.write(node.as_csv() + ',,F2\n')
    out_f.close()
    out_f = open(os.path.join(output_dir, "edges.csv"), "w")
    for edge in node_holder.get_edges():
        out_f.write(edge.as_csv() + ',FRAG\n')
    out_f.close()


def build_network(attrs, node_holder, max_frags=0, base_dir='.', verbosity=0, recurse=True):

    log_file = None
    if ENABLE_BUILD_NETWORK_LOG:
        log_file_name = os.path.join(base_dir, BUILD_NETWORK_LOG)
        log_file = open(log_file_name, 'w')

    # Create the nodes and test with output
    tqdm_disable = True #if verbosity else False
    for attr in tqdm(attrs, disable=tqdm_disable):

        start_time = timeit.default_timer()
        node, is_new = node_holder.create_or_retrieve_node(attr.SMILES)
        retrieve_end_time = timeit.default_timer()
        create_end_time = None
        direct_frags = 0
        if node and is_new:
            direct_frags, direct_ring_ring_splits =\
                create_children(node, node_holder,
                                max_frags, attr.SMILES,
                                log_file, recurse=recurse)
            #print("Direct frags {0}".format(str(node_holder.size())))
            create_end_time = timeit.default_timer()

        if verbosity:
            total_dur = create_end_time - start_time
            print('{} {} {}'.format(attr.SMILES, total_dur, direct_ring_ring_splits + direct_frags))

        if log_file:
            nh_nodes, nh_edges = node_holder.size()
            node_created = True if node else False
            if is_new:
                retrieve_dur = retrieve_end_time - start_time
                create_dur = create_end_time - retrieve_end_time
                total_dur = create_end_time - start_time
                log_file.write('1,%s,%s,%s,%s,%s,%s,%s,%s\n' % (attr.SMILES, node_created, retrieve_dur, create_dur, total_dur, nh_nodes, nh_edges, direct_frags))
                log_file.flush()
                os.fsync(log_file)
            else:
                retrieve_dur = retrieve_end_time - start_time
                log_file.write('0,%s,%s,%s,%s,%s,%s,%s\n' % (attr.SMILES, node_created, retrieve_dur, 0, retrieve_dur, nh_nodes, nh_edges))
                log_file.flush()
                os.fsync(log_file)

    if log_file:
        log_file.close()

    return node_holder


def add_node(tx, smiles, hac, chac, osmiles):
    tx.run(
        "MERGE (:F2 { smiles: $smiles, hac: toInt($hac), chac: toInt($chac), osmiles: $osmiles})",
        smiles=smiles,
        hac=hac,
        chac=chac,
        osmiles=osmiles,
    )


def add_edge(tx, smiles, smiles_two, edge_meta):
    tx.run(
        "MATCH (n1:F2 { smiles: $smiles}), (n2:F2 { smiles: $smiles_two}) MERGE (n1)-[:FRAG{label:$edge_meta}]->(n2)",
        smiles=smiles,
        smiles_two=smiles_two,
        edge_meta=edge_meta,
    )


def add_attr(tx, smiles, attr):
    tx.run(
        "MATCH (n:F2 { smiles: $smiles} ) set n:MOL, n:Mol, n.Mol=$attr",
        smiles=smiles,
        attr=attr,
    )
