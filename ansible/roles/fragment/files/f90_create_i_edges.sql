/*
 * Create Fragmentation Database SQL Statements: 
 * Purpose: Creates fragmentation company specific tables
 *
 * Author | Date    | Version
 * Duncan | 03/2020 | Initial Version
 *
 */

/*
 * Create i_edge
 */

DROP TABLE IF EXISTS i_edge;
CREATE TABLE i_edge (
  p_smiles TEXT,
  c_smiles TEXT,
  present BOOLEAN,
  label TEXT
);

