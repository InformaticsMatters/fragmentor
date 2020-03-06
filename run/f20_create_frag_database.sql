/*
 * Create Fragmentation Database SQL Statements: 
 * Purpose: Creates fragmentation company specific tables
 *
 * Called from p10_create_ISO_database.sh
 *
 * Author | Date    | Version
 * Duncan | 03/2020 | Initial Version
 *
 */


/*
 * Drop all tables if they exist (NB Order because of table key constrants) 
 */

DROP TABLE IF EXISTS i_edge;

DROP TABLE IF EXISTS i_mols;

DROP TABLE IF EXISTS i_node;

/*
 * Create i_edge
 */
CREATE TABLE i_edge (
  p_smiles TEXT,
  c_smiles TEXT,
  present BOOLEAN,
  label TEXT
);

/* 
 * Create i_mols        
 */
CREATE TABLE i_mols (
  osmiles TEXT,
  isosmiles TEXT,
  nonisosmiles TEXT,
  hac SMALLINT,
  rac SMALLINT,
  cmpd_id TEXT,
  isomol_id INTEGER,
  nonisomol_id INTEGER
);


/* 
 * Create i_node             
 */
CREATE TABLE i_node (
  smiles TEXT,
  nonisomol_id INTEGER,
  hac SMALLINT,
  rac SMALLINT,
  child_count INTEGER,
  edge_count INTEGER,
  ftime INTEGER
);
