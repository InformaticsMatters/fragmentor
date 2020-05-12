/*
 * Extract Neo4j Edge Data SQL Statements:
 * Purpose: Extract Neo4j Edge data into CSV file on database server
 */

--
-- Edges.csv
-- :START_ID(F2),:END_ID(F2),label,:TYPE
--

COPY (select distinct np.smiles, nc.smiles, label, 'FRAG'
        from o_edge_parent os
        join pg_class p on os.tableoid = p.oid and p.relname in (%(SOURCETABLES)s)
        join nonisomol np ON np.id = os.parent_id
        join nonisomol nc ON nc.id = os.child_id)
  TO %(NEOEDGEFILE)s DELIMITER ',' CSV;

--SELECT distinct c.parent_id, c.child_id, c.label
--  FROM o_edge_parent c, pg_class p
-- WHERE c.tableoid = p.oid and p.relname in ('o_edge_xchem_dsip','o_edge_molport');
