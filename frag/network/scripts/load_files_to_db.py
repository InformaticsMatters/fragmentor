"""
Insert data from nodes and edges files into the Neo4j database.
The logic is as follows:

Read the nodes.csv file and for each line:
- check whether the SMILES (first column) already has a node in the db
  - if so then no further action
  - if not then
    - insert a new node with the data for that molecule
    - add the SMILES to the set that have been inserted
- if any nodes have been inserted then process the edges.csv:
  - for each line:
    - if the 'to' or 'from' SMILES is present in the set of inserted SMILES then this will be a new edge so:
      - insert the edge with the corresponding data
    - if not then this edge between the 2 SMILES was already in the database


Typical node line:
CC(C)(C)OC(=O)N1CCCC1,12,5,CC(C)(C)OC(O)C1CCCC1,LPQZERIRKRYGGM-UHFFFAOYNA-N,"InChI=1/C9H17NO2/c1-9(2,3)12-8(11)10-6-4-5-7-10/h4-7H2,1-3H3MA",,F2

Typical edge line:
CC(C)(C)OC(=O)N1CCCC1,C1CCNC1,FG|CC(C)(C)OC(=O)[Xe]|CC(C)(C)OC(O)[100Xe]|RING|[Xe]N1CCCC1|[100Xe]C1CCCC1,FRAG

TODO:
- add the isomol data
- add the supplier data
- patch in the cmpd_id property for existing nodes
"""

import argparse, os, time
import csv

from neo4j import GraphDatabase


def load_files(nodes, edges, supplier, neo4j_server, neo4j_username, neo4j_password, rollback=False, load_id=None):
    if not os.path.isfile(nodes):
        print(nodes, 'does not exist')
        exit(1)
    if not os.path.isfile(edges):
        print(edges, 'does not exist')
        exit(1)

    print('Processing', nodes, edges, supplier)

    t0 = time.time()

    driver = GraphDatabase.driver(neo4j_server, auth=(neo4j_username, neo4j_password))
    with driver.session() as session:
        with session.begin_transaction() as tx:
            new_nodes = load_nodes(tx, nodes, load_id)
            if new_nodes:
                new_edge_count = load_edges(tx, edges, new_nodes, load_id)

            if rollback:
                session.close()
                print('Transaction rolled back. No changes made.')

    t1 = time.time()
    print('Took', (t1 - t0), 'seconds')


def load_nodes(tx, nodes, load_id):
    count = 0
    hits = 0
    misses = 0
    loaded = set()
    with open(nodes, 'rt') as csvfile:
        reader = csv.reader(csvfile)
        for tokens in reader:
            count += 1
            # smiles:ID(F2),hac:INT,chac:INT,osmiles,inchik,inchis,cmpd_ids:STRING[],:LABEL
            smiles = tokens[0]
            hac = int(tokens[1])
            rac = int(tokens[2])
            osmiles = tokens[3]
            inchik = tokens[4]
            inchis = tokens[5]
            cmpd_ids = tokens[6]
            labels = tokens[7]

            found, inserted = process_node(tx, smiles, hac, rac, osmiles, inchik, inchis, cmpd_ids, labels, load_id)
            if found:
                hits += 1
            else:
                misses += 1
                loaded.add(smiles)
                if not inserted:
                    print('WARNING: node', smiles, 'not inserted')

    print('Processed', count, 'nodes', hits, 'found', misses, 'missing')
    return loaded


def process_node(tx, smiles, hac, rac, osmiles, inchik, inchis, cmpd_ids, labels, load_id):
    query = 'MATCH (n:F2 { smiles: $smi }) return n.smiles'
    result = tx.run(query, smi=smiles)
    inserted = False
    if result.single():
        found = True
    else:
        found = False
        insert = 'CREATE (n:' + labels.replace(';', ':') + \
                 ' { smiles: $smi, osmiles: $osmi, inchik: $inchik, inchis: $inchis, hac: $hac, chac: $rac, cmpd_ids: $ids'
        if load_id:
            insert = insert + ', load_id: $load_id'
        insert = insert + ' }) RETURN n'

        # unclear whether cmpd_ids can have multiple value
        ids = [cmpd_ids]
        result = tx.run(insert, smi=smiles, osmi=osmiles, inchik=inchik, inchis=inchis, hac=hac, rac=rac, ids=ids, load_id=load_id)
        record = result.single()
        if record:
            inserted = True

    return found, inserted


def load_edges(tx, edges, new_nodes, load_id):
    count = 0
    insert_count = 0
    skipped_count = 0
    failed_count = 0
    with open(edges, 'rt') as csvfile:
        reader = csv.reader(csvfile)
        for tokens in reader:
            count += 1
            n1 = tokens[0]
            n2 = tokens[1]
            label = tokens[2]
            type = tokens[3]
            if type != 'FRAG':
                print('Unexpected relationship type:', type)
                exit(1)
            if n1 in new_nodes or n2 in new_nodes:
                stmt = 'MATCH (n1:F2 { smiles: $n1 }), (n2:F2 { smiles: $n2 })' + \
                       ' CREATE (n1)-[r:FRAG { label: $label'
                if load_id:
                    stmt = stmt + ', load_id: $load_id'
                stmt = stmt + ' }]->(n2) RETURN r'
                result = tx.run(stmt, n1=n1, n2=n2, label=label, load_id=load_id)
                if result.single():
                    insert_count += 1
                else:
                    failed_count += 1
            else:
                skipped_count += 1
    print('Inserted', insert_count, 'edges,', skipped_count, 'skipped,', failed_count, 'failed')
    return insert_count


def main():

    # Example (combine):

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Load nodes and edges files')
    parser.add_argument('-i', '--input-dir', required=True, help="Directory where inputs are found")
    parser.add_argument('-s', '--supplier-name', required=True, help="Name of the supplier")
    parser.add_argument('-r', '--rollback', action='store_true', help="Rollback the transaction so that no changes are made")
    parser.add_argument('-l', '--load-id', help="Add this value as the 'load_id' property to the inserted nodes and edges")
    parser.add_argument('--server', help='Neo4j database server')
    parser.add_argument('--username', help='Neo4j username')
    parser.add_argument('--password', help='Neo4j password')

    args = parser.parse_args()
    print("load_files_to_db: ", args)

    dir = args.input_dir
    if not os.path.isdir(dir):
        print("--input-dir argument is not a directory")
        exit(1)

    if args.server:
        neo4j_server = args.server
    else:
        neo4j_server = os.getenv('NEO4J_SERVER')
    if not neo4j_server:
        print('WARNING: no database server URL provided')
        exit(1)

    if args.username:
        neo4j_username = args.username
    else:
        neo4j_username = os.getenv('NEO4J_USERNAME')
    if not neo4j_username:
        print('WARNING: no database username provided')
        exit(1)

    if args.password:
        neo4j_password = args.password
    else:
        neo4j_password = os.getenv('NEO4J_PASSWORD')
    if not neo4j_password:
        print('WARNING: no database password provided')
        exit(1)

    load_files(os.path.join(dir, 'nodes.csv'), os.path.join(dir, 'edges.csv'), args.supplier_name,
               neo4j_server, neo4j_username, neo4j_password,
               rollback=args.rollback, load_id=args.load_id)


if __name__ == "__main__":
    main()
