/*
 * Load Fragmented Nodes Results SQL Statements: 
 * Purpose: Loads fragmented node data into Frag database
 *
 * Called from p80_load_frag_results.sh
 *
 * Author | Date    | Version
 * Duncan | 03/2020 | Initial Version
 *
 */

/*
 * Update i_node with smailes from existing nonisomol table 
 */
\timing

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

INSERT INTO nonisomol (smiles, hac, rac, child_count, edge_count, ftime)
  SELECT smiles, hac, rac, child_count, edge_count, ftime FROM i_node
  WHERE nonisomol_id IS NULL;
commit;

/*
 * Update counts 
 */
begin;
SELECT clock_timestamp();

-- this next query needs optimising
WITH s AS (SELECT i.nonisomol_id, i.child_count, i.edge_count, i.ftime FROM i_node i
  JOIN nonisomol n ON n.id = i.nonisomol_id
  WHERE i.nonisomol_id = n.id AND i.nonisomol_id IS NOT NULL)
UPDATE nonisomol n SET child_count = s.child_count, edge_count = s.edge_count, ftime = s.ftime FROM s
  WHERE s.nonisomol_id = n.id;
commit;
