/*
 * Extract Neo4j Node Data SQL Statements:
 * Purpose: Extract Neo4j node data into CSV file on database server
 */

--
-- Nodes.csv
-- smiles:ID(F2),hac:INT,chac:INT,osmiles,cmpd_ids:STRING[],:LABEL
--

COPY (WITH RECURSIVE fragments AS (
        select parent_id, child_id, parent_smiles, child_smiles, hac, rac, ring_smiles, child_hac, child_rac, child_ring_smiles, child_child_count
          from v_edge_node e
         inner join mol_source ms on e.parent_id = ms.nonisomol_id and ms.source_id = %(SOURCEID)s
        union
         select c.parent_id, c.child_id, c.parent_smiles, c.child_smiles, c.hac, c.rac, c.ring_smiles, c.child_hac, c.child_rac, c.child_ring_smiles,c.child_child_count
           from v_edge_node c
          inner join fragments p on c.parent_id = p.child_id
        ) select f1.parent_smiles, f1.hac, f1.rac, f1.ring_smiles, ms.code as cmpd,
            (case when ms.code notnull then 'F2;CanSmi;Mol;V_XCHEM'
                  else 'F2'
             end) as label
         from fragments f1
         left join mol_source ms on f1.parent_id = ms.nonisomol_id and ms.source_id = %(SOURCEID)s
        union
          select f2.child_smiles, f2.child_hac, f2.child_rac, f2.child_ring_smiles, NULL, 'F2'
         from fragments f2
        where  f2.child_child_count = 0)
         TO %(NEONODEFILE)s DELIMITER ',' CSV;
