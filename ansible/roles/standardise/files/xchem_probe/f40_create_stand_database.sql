-- Create Fragmentation Database SQL Statements for Libary xchem_probe
-- Purpose: Creates fragmentation company specific tables. This is the same as the xchem_dsip library.

DROP TABLE IF EXISTS i_mols_dsip;

-- Create i_mols
CREATE TABLE i_mols_dsip (
  osmiles TEXT,
  isosmiles TEXT,
  nonisosmiles TEXT,
  hac SMALLINT,
  cmpd_id TEXT,
  isomol_id INTEGER,
  nonisomol_id INTEGER
);
