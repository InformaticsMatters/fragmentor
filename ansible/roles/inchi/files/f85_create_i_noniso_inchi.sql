/*
 * Create Fragmentation i_noniso_inchi Table SQL Statements :
 * Purpose: Creates table to load inchi data for new nonisomol records
 *
 * Called from f85_load_noniso_inchi.sh
 *
 * Author | Date    | Version
 * Duncan | 03/2020 | Initial Version
 *
 */

/*
 * Create i_node             
 */
DROP TABLE IF EXISTS i_noniso_inchi;
CREATE TABLE i_noniso_inchi (
  smiles TEXT,
  sinchik TEXT,
  sinchis TEXT,
  ninchik TEXT,
  ninchis TEXT
);