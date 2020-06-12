-- Purpose: Loads chemspace standardised data into database

-- Load nonisomol except where already exists
begin;
INSERT INTO nonisomol (smiles, hac, source_id)
  SELECT nonisosmiles, hac, %(SOURCEID)s from i_mols_chemspace
  ON CONFLICT ON CONSTRAINT nonisomol_smiles_key DO NOTHING;
commit;

-- Update i_mols_chemspace with nonisomol_id
begin;
UPDATE i_mols_chemspace i SET nonisomol_id = n.id
  FROM nonisomol n
    WHERE n.smiles = i.nonisosmiles;
commit;

-- Load isomol from i_mols_chemspace
begin;
INSERT INTO isomol (smiles, nonisomol_id, source_id)
  SELECT isosmiles, nonisomol_id, %(SOURCEID)s from i_mols_chemspace
  WHERE isosmiles != nonisosmiles
  ON CONFLICT ON CONSTRAINT isomol_smiles_key DO NOTHING;
commit;

-- Update i_mols with isomol_id
begin;
UPDATE i_mols_chemspace i SET isomol_id = iso.id
  FROM isomol iso
    WHERE iso.smiles = i.isosmiles;
commit;

-- Load mol_source from i_mols
begin;
INSERT INTO mol_source (smiles, code, source_id, lead_time, nonisomol_id, isomol_id)
  (SELECT osmiles, cmpd_id, %(SOURCEID)s, 0, nonisomol_id, isomol_id FROM i_mols_chemspace);
commit;

-- Load price information from i_mols and mol_source
begin;
INSERT INTO price (quantity_mg, price, price_min, price_max, molsource_id)
  (SELECT 1, im.price,  im.price,  im.price, ms.id
   FROM  i_mols_chemspace im, mol_source ms
   WHERE im.nonisomol_id = ms.nonisomol_id
     AND ms.source_id =%(SOURCEID)s
  );
commit;


