/*
 * Create Fragmentation Database SQL Statements: 
 * Purpose: Creates fragmentation company specific tables
 *
 * Called from p80_load_frag_results.sh
 *
 * Author | Date    | Version
 * Duncan | 03/2020 | Initial Version
 *
 */

/*
 * Create i_edge
 */

DROP INDEX IF EXISTS ix_i_edge_psmiles;
DROP INDEX IF EXISTS ix_i_edge_csmiles;

DROP TABLE IF EXISTS i_edge;
CREATE TABLE i_edge (
  p_smiles TEXT,
  c_smiles TEXT,
  present BOOLEAN,
  label TEXT
);

