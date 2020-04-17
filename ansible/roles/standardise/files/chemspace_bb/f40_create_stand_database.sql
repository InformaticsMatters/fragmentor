/*
 * Create Fragmentation Database SQL Statements for Vendor CHEMSPACE
 * Purpose: Creates fragmentation company specific tables
 *
 * Author | Date    | Version
 * Duncan | 03/2020 | Initial Version
 *
 */


/*
 * Drop all tables if they exist (NB Order because of table key constrants) 
 */

DROP TABLE IF EXISTS i_mols_chemspace;

/*
 * Create i_mols_chemspace - has price, but no rac (in comparison to dsip)
 */
CREATE TABLE i_mols_chemspace (
  osmiles TEXT,
  isosmiles TEXT,
  nonisosmiles TEXT,
  hac SMALLINT,
  cmpd_id TEXT,
  price INTEGER,
  isomol_id INTEGER,
  nonisomol_id INTEGER
);
