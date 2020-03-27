/*
 * Create Fragmentation i_node Table SQL Statements :
 * Purpose: Creates fragmentation company specific tables
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
  child_count INTEGER,
  edge_count INTEGER,
  ftime INTEGER
);
