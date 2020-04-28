/*
 * Extract Neo4j Edge Data SQL Statements:
 * Purpose: Extract Neo4j Edge data into CSV file on database server
 */

COPY (with RECURSIVE fragments AS (
       select parent_smiles, child_smiles, label, 'FRAG', parent_id, child_id
         from v_edge e
         inner join mol_source ms on e.parent_id = ms.nonisomol_id and ms.source_id = %(SOURCEID)s
         union
            select c.parent_smiles, c.child_smiles, c.label, 'FRAG', c.parent_id, c.child_id
             from v_edge c
             inner join fragments p on c.parent_id = p.child_id
        ) select parent_smiles, child_smiles, label, 'FRAG' from fragments)
        TO %(OUTFILE)s DELIMITER ',' CSV;
