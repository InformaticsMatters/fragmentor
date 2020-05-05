/*
 * Load Fragmented Nodes Results SQL Statements: 
 * Purpose: Loads fragmented node data into Frag database
 *
 * Author | Date    | Version
 * Duncan | 03/2020 | Initial Version
 *
 */

/*
 * Update isomol from i_node with the inchi key.
 */
begin;
SELECT clock_timestamp();

WITH s
    AS (SELECT i.smiles, i.ninchik, i.ninchis
  FROM i_iso_inchi i
  JOIN isomol iso ON iso.smiles = i.smiles
 WHERE i.ninchik <> '')
UPDATE isomol iso
   SET inchik = s.ninchik, inchis = s.ninchis
  FROM s
 WHERE s.smiles = iso.smiles;
commit;
