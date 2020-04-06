/*
 * Create Fragmentation i_inchi Table SQL Statements :
 * Purpose: Creates table to load inchi data
 *
 * Called from f85_load_inchi.sh
 *
 * Author | Date    | Version
 * Duncan | 03/2020 | Initial Version
 *
 */



/*
 * Create i_node             
 */
DROP TABLE IF EXISTS i_inchi;
CREATE TABLE i_inchi (
  smiles TEXT,
  sinchik TEXT,
  sinchis TEXT,
  ninchik TEXT,
  ninchis TEXT
);
