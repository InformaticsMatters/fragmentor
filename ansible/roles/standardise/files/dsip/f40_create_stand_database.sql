/*
 * Create Fragmentation Database SQL Statements for Vendor DSIP
 * Purpose: Creates fragmentation company specific tables
 *
 * Called from f40_create_stand_database.sh
 *
 * Author | Date    | Version
 * Duncan | 03/2020 | Initial Version
 *
 */


/*
 * Drop all tables if they exist (NB Order because of table key constrants) 
 */

DROP TABLE IF EXISTS i_mols_dsip;

/*
 * Create i_mols        
 */
CREATE TABLE i_mols_dsip (
  osmiles TEXT,
  isosmiles TEXT,
  nonisosmiles TEXT,
  hac SMALLINT,
  cmpd_id TEXT,
  isomol_id INTEGER,
  nonisomol_id INTEGER
);
