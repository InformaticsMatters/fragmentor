/*
 * Create Fragmentation i_node Table SQL Statements :
 * Purpose: Creates fragmentation specific i_node table
 *
 * Called from p80_load_nodes.sh
 *
 * Author | Date    | Version
 * Duncan | 03/2020 | Initial Version
 *
 */



/*
 * Create i_node             
 */
DROP TABLE IF EXISTS i_node;
CREATE TABLE i_node (
  smiles TEXT,
  nonisomol_id INTEGER,
  hac SMALLINT,
  rac SMALLINT,
  ring_smiles TEXT,
  child_count INTEGER,
  edge_count INTEGER,
  ftime INTEGER
);
