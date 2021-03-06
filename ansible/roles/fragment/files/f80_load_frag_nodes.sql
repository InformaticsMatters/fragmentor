/*
 * Load Fragmented Nodes Results SQL Statements: 
 * Purpose: Loads fragmented node data into Frag database
 *
 * Author | Date    | Version
 * Duncan | 03/2020 | Initial Version
 *
 */

begin;
SELECT clock_timestamp();

UPDATE i_node i SET nonisomol_id = n.id
  FROM nonisomol n WHERE n.smiles = i.smiles;
commit;

/*
 * Insert new Nodes into nonismol
 */
begin;
SELECT clock_timestamp();

INSERT INTO nonisomol (smiles, hac, rac, ring_smiles, child_count, edge_count, ftime, source_id)
  SELECT smiles, hac, rac, ring_smiles, child_count, edge_count, ftime, %(SOURCEID)s FROM i_node
  WHERE nonisomol_id IS NULL;
commit;

/*
 * Update counts
 */
begin;
SELECT clock_timestamp();

-- this next query needs optimising
WITH s AS (SELECT i.nonisomol_id, i.hac, i.rac, i.ring_smiles, i.child_count, i.edge_count, i.ftime FROM i_node i
  JOIN nonisomol n ON n.id = i.nonisomol_id
  WHERE i.nonisomol_id = n.id AND i.nonisomol_id IS NOT NULL)
UPDATE nonisomol n SET hac = s.hac, rac = s.rac, ring_smiles = s.ring_smiles, child_count = s.child_count, edge_count = s.edge_count, ftime = s.ftime FROM s
  WHERE s.nonisomol_id = n.id;
commit;

