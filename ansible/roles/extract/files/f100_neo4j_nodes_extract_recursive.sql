/*
 * Extract Neo4j Node Data SQL Statements:
 * Purpose: Extract Neo4j node data into CSV file on database server
 */

--
-- Nodes.csv
-- smiles:ID(F2),hac:INT,chac:INT,osmiles,cmpd_ids:STRING[],:LABEL
--

COPY (WITH RECURSIVE fragments AS (
           select parent_id, child_id, label
             from edge e
            inner join mol_source ms on e.parent_id = ms.nonisomol_id and ms.source_id = %(SOURCEID)s
            union
           select c.parent_id, c.child_id, c.label
             from edge c
             inner join fragments p on c.parent_id = p.child_id
        ) select np.smiles, np.hac, np.rac, np.ring_smiles, np.inchik, np.inchis, ms.code as cmpd,
           (case when ms.code notnull then 'F2;CanSmi;Mol;' || v.supplier_node_name
              else 'F2'
            end) as label
            from fragments os
            join nonisomol np ON np.id = os.parent_id
            left join mol_source ms on os.parent_id = ms.nonisomol_id and ms.source_id = %(SOURCEID)s
            left join source s on s.id = ms.source_id
            left join vendor_name v on s.name = v.vendor_name
           union
          select nc.smiles, nc.hac, nc.rac, nc.ring_smiles, nc.inchik, nc.inchis, NULL, 'F2'
            from fragments os2
            join nonisomol nc ON nc.id = os2.child_id
             and nc.child_count = 0
     ) TO %(NEONODEFILE)s DELIMITER ',' CSV;

