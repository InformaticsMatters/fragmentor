/*
 * Load Standardised Data SQL Statements: 
 * Purpose: Loads standardised data into ISO database
 *
 * Called from p40_load_standardised_data.sh
 */

/*
 * Load nonisomol except where already exists 
 */
begin;
INSERT INTO nonisomol (smiles, hac, rac)
  SELECT nonisosmiles, hac, rac from i_mols
  ON CONFLICT ON CONSTRAINT nonisomol_smiles_key DO NOTHING;
commit;

/*
 * Update i_mols with nonisomol_id 
 */
begin;
UPDATE i_mols i SET nonisomol_id = n.id
  FROM nonisomol n
    WHERE n.smiles = i.nonisosmiles;
commit;

/*
 * Load isomol from i_mols 
 */
begin;
INSERT INTO isomol (smiles, nonisomol_id)
  SELECT isosmiles, nonisomol_id from i_mols
  WHERE isosmiles != nonisosmiles
  ON CONFLICT ON CONSTRAINT isomol_smiles_key DO NOTHING;
commit;

/*
 * Update i_mols with isomol_id 
 */
begin;
UPDATE i_mols i SET isomol_id = iso.id
  FROM isomol iso
    WHERE iso.smiles = i.isosmiles;
commit;

/*
 * Load mol_source from i_mols 
 */
begin;
INSERT INTO mol_source (smiles, code, source_id, nonisomol_id, isomol_id)
  (SELECT osmiles, cmpd_id, 6, nonisomol_id, isomol_id FROM i_mols);
commit;
