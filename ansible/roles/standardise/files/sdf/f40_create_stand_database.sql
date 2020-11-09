/*
 * Create Fragmentation Database SQL Statements for Vendor SDF
 * Purpose: Creates fragmentation company specific tables
 *
 * Author | Date    | Version
 * Duncan | 03/2020 | Initial Version
 *
 */


/*
 * Drop all tables if they exist (NB Order because of table key constrants) 
 */

DROP TABLE IF EXISTS i_mols_sdf;

/*
 * Create i_sdf
 */
CREATE TABLE i_mols_sdf (
  osmiles TEXT,
  isosmiles TEXT,
  nonisosmiles TEXT,
  hac SMALLINT,
  cmpd_id TEXT,
  isomol_id INTEGER,
  nonisomol_id INTEGER
);
