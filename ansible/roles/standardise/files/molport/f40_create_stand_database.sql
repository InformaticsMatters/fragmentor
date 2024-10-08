/*
 * Create Fragmentation Database SQL Statements for Vendor Molport
 * Purpose: Creates fragmentation company specific tables
 *
 * Author | Date    | Version
 * Alan   | 10/2024 | BIGINT Version
 * Duncan | 03/2020 | Initial Version
 *
 */

DROP TABLE IF EXISTS i_mols_molport;

/*
 * Create i_mols_molport - has price columns, no rac and lead time (in comparison to dsip)
 */
CREATE TABLE i_mols_molport (
  osmiles TEXT,
  isosmiles TEXT,
  nonisosmiles TEXT,
  hac SMALLINT,
  cmpd_id TEXT,
  price_1 TEXT,
  price_5 TEXT,
  price_50 TEXT,
  lead_time TEXT,
  isomol_id bigint,
  nonisomol_id bigint
);
