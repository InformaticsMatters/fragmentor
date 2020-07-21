-- Purpose: Loads molport standardised data into database

-- Load nonisomol except where already exists
begin;
INSERT INTO nonisomol (smiles, hac, source_id)
  SELECT nonisosmiles, hac, %(SOURCEID)s from i_mols_molport
  ON CONFLICT ON CONSTRAINT nonisomol_smiles_key DO NOTHING;
commit;

-- Update i_mols with nonisomol_id
begin;
UPDATE i_mols_molport i SET nonisomol_id = n.id
  FROM nonisomol n
    WHERE n.smiles = i.nonisosmiles;
commit;

-- Load isomol from i_mols
begin;
INSERT INTO isomol (smiles, nonisomol_id, source_id)
  SELECT isosmiles, nonisomol_id, %(SOURCEID)s from i_mols_molport
  WHERE isosmiles != nonisosmiles
  ON CONFLICT ON CONSTRAINT isomol_smiles_key DO NOTHING;
commit;

-- Update i_mols with isomol_id
begin;
UPDATE i_mols_molport i SET isomol_id = iso.id
  FROM isomol iso
    WHERE iso.smiles = i.isosmiles;
commit;

-- Load mol_source from i_mols
-- Note: <> 'null' - it would be better to use force_null on the copy command, but this may affect other vendors.
begin;
INSERT INTO mol_source (smiles, code, source_id, lead_time, nonisomol_id, isomol_id)
  (SELECT osmiles, cmpd_id, %(SOURCEID)s,
         case
            when lead_time <> 'null' then CAST (lead_time AS INT)
            else 0
         end,
         nonisomol_id, isomol_id
  FROM i_mols_molport);
commit;

-- Load price information from i_mols and mol_source
begin;
INSERT INTO price (quantity_mg, price, price_min, price_max, molsource_id)
  (
   SELECT 1, 0,
       case
          when SUBSTRING(im.price_1, '- ') notnull then CAST(split_part(im.price_1, '-', 1) as INT)
          else 0
       end,
       case
          when SUBSTRING(im.price_1, '< ') notnull then CAST(split_part(im.price_1, '<', 2) as INT)
          when SUBSTRING(im.price_1, '- ') notnull then CAST(split_part(im.price_1, '-', 2) as INT)
          else 0
       end,
       ms.id
   from  i_mols_molport im, mol_source ms
   where (im.nonisomol_id = ms.nonisomol_id and im.isomol_id = ms.isomol_id)
       or (im.nonisomol_id = ms.nonisomol_id and im.isomol_id is null and ms.isomol_id is null)
     and ms.source_id = %(SOURCEID)s
   union all
   SELECT 5, 0,
       case
          when SUBSTRING(im.price_5, '- ') notnull then CAST(split_part(im.price_5, '-', 1) as INT)
          else 0
       end,
       case
          when SUBSTRING(im.price_5, '< ') notnull then CAST(split_part(im.price_5, '<', 2) as INT)
          when SUBSTRING(im.price_5, '- ') notnull then CAST(split_part(im.price_5, '-', 2) as INT)
          else 0
       end,
       ms.id
   from  i_mols_molport im, mol_source ms
   where (im.nonisomol_id = ms.nonisomol_id and im.isomol_id = ms.isomol_id)
       or (im.nonisomol_id = ms.nonisomol_id and im.isomol_id is null and ms.isomol_id is null)
     and ms.source_id = %(SOURCEID)s
   union all
   SELECT 50, 0,
       case
          when SUBSTRING(im.price_50, '- ') notnull then CAST(split_part(im.price_50, '-', 1) as INT)
          else 0
       end,
       case
          when SUBSTRING(im.price_50, '< ') notnull then CAST(split_part(im.price_50, '<', 2) as INT)
          when SUBSTRING(im.price_50, '- ') notnull then CAST(split_part(im.price_50, '-', 2) as INT)
          else 0
       end,
       ms.id
   from  i_mols_molport im, mol_source ms
   where (im.nonisomol_id = ms.nonisomol_id and im.isomol_id = ms.isomol_id)
       or (im.nonisomol_id = ms.nonisomol_id and im.isomol_id is null and ms.isomol_id is null)
     and ms.source_id = %(SOURCEID)s
  );
commit;
