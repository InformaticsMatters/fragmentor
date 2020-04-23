/*
 * Load Standardised Data SQL Statements: 
 * Purpose: Loads enamine standardised data into frag database
 */

/*
 * Load nonisomol except where already exists 
 */
begin;
SELECT clock_timestamp();

INSERT INTO nonisomol (smiles, hac, rac, source_id)
  SELECT nonisosmiles, hac, rac, %(SOURCEID)s from i_mols_enamine
  ON CONFLICT ON CONSTRAINT nonisomol_smiles_key DO NOTHING;
commit;

/*
 * Update i_mols with nonisomol_id 
 */

begin;
SELECT clock_timestamp();

UPDATE i_mols_enamine i SET nonisomol_id = n.id
  FROM nonisomol n
    WHERE n.smiles = i.nonisosmiles;
commit;

/*
 * Load isomol from i_mols 
 */

begin;
SELECT clock_timestamp();

INSERT INTO isomol (smiles, nonisomol_id, source_id)
  SELECT isosmiles, nonisomol_id, %(SOURCEID)s from i_mols_enamine
  WHERE isosmiles != nonisosmiles
  ON CONFLICT ON CONSTRAINT isomol_smiles_key DO NOTHING;
commit;

/*
 * Update i_mols with isomol_id 
 */

begin;
SELECT clock_timestamp();

UPDATE i_mols_enamine i SET isomol_id = iso.id
  FROM isomol iso
    WHERE iso.smiles = i.isosmiles;
commit;

/*
 * Load mol_source from i_mols 
 */

begin;
SELECT clock_timestamp();

INSERT INTO mol_source (smiles, code, source_id, lead_time, nonisomol_id, isomol_id)
  (SELECT osmiles, cmpd_id, %(SOURCEID)s, 0, nonisomol_id, isomol_id FROM i_mols_enamine);
commit;
