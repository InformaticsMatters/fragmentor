/*
 * Create Fragmentation Database SQL Statements for Vendor Enamine
 * Purpose: Creates fragmentation company specific tables
 *
 * Author | Date    | Version
 * Alan   | 10/2024 | BIGINT Version
 * Duncan | 03/2020 | Initial Version
 *
 */


/*
 * Drop all tables if they exist (NB Order because of table key constrants)
 */

DROP TABLE IF EXISTS i_mols_enamine;

/*
 * Create i_mols
 */
CREATE TABLE i_mols_enamine (
  osmiles TEXT,
  isosmiles TEXT,
  nonisosmiles TEXT,
  hac SMALLINT,
  rac SMALLINT,
  cmpd_id TEXT,
  isomol_id bigint,
  nonisomol_id bigint
);
