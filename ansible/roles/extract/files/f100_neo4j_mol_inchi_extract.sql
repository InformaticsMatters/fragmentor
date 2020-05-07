/*
 * Extract Neo4j Molecule Inchi Edges Data SQL Statements:
 * Purpose: Extract Neo4j Molecule Inchi Edges data into CSV file on database server
 */

--
-- molecule-inchi-edges.csv
-- :START_ID(NONISOMOL),:END_ID(Inchi),:TYPE
--- Link inchi-nodes to nonisomol-nodes, "HasInchi"
--
COPY (select np.smiles, np.inchi_id, 'HasInchi' as label
       from o_source_edge os
       join nonisomol np ON np.id = os.parent_id
       left join mol_source ms on os.parent_id = ms.nonisomol_id and ms.source_id = %(SOURCEID)s
      union
     select nc.smiles, nc.inchi_id, 'HasInchi' as label
       from o_source_edge os2
       join nonisomol nc ON nc.id = os2.child_id
        and nc.child_count = 0)
   TO %(MOLINCHIEDGEFILE)s DELIMITER ',' CSV;
