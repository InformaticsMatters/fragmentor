/*
 * Extract Neo4j Edge Data SQL Statements:
 * Purpose: Extract Neo4j Edge data into CSV file on database server
 */

--
-- Edges.csv
-- :START_ID(F2),:END_ID(F2),label,:TYPE
--

COPY (select np.smiles, nc.smiles, label, 'FRAG'
        from o_source_edge os
        join nonisomol np ON np.id = os.parent_id
        join nonisomol nc ON nc.id = os.child_id)
  TO %(NEOEDGEFILE)s DELIMITER ',' CSV;

