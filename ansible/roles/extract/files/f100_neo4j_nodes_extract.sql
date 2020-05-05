/*
 * Extract Neo4j Node Data SQL Statements:
 * Purpose: Extract Neo4j node data into CSV file on database server
 */

--
-- Nodes.csv
-- smiles:ID(F2),hac:INT,chac:INT,osmiles,cmpd_ids:STRING[],:LABEL
--
COPY (select np.smiles, np.hac, np.rac, np.ring_smiles, ms.code as cmpd,
        (case when ms.code notnull then 'F2;CanSmi;Mol;V_XCHEM'
              else 'F2'
         end) as label
       from o_source_edge os
       join nonisomol np ON np.id = os.parent_id
       left join mol_source ms on os.parent_id = ms.nonisomol_id and ms.source_id = 4
      union
     select nc.smiles, nc.hac, nc.rac, nc.ring_smiles, NULL, 'F2'
      from o_source_edge os2
       join nonisomol nc ON nc.id = os2.child_id
        and nc.child_count = 0)
   TO %(NEONODEFILE)s DELIMITER ',' CSV;
