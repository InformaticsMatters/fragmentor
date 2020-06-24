/*
 * Load Nonisomol Inchi Results SQL Statements:
 * Purpose: Loads Inchi key data into Frag database
 *
 * Author | Date    | Version
 * Duncan | 03/2020 | Initial Version
 *
 */

/*
 * Insert INCHI tables.
 */
begin;
INSERT INTO inchi (inchik, inchis)
  SELECT sinchik, sinchis from i_noniso_inchi
  ON CONFLICT ON CONSTRAINT inchi_inchik_key DO NOTHING;
commit;

/*
 * Update inchi key for nonisomols.
 */
begin;
WITH s
    AS (SELECT i.smiles, i.ninchik, i.ninchis, inc.id
  FROM i_noniso_inchi i
  JOIN inchi inc ON inc.inchik = i.sinchik
 WHERE i.ninchik <> '')
UPDATE nonisomol non
   SET inchik = s.ninchik, inchis = s.ninchis, inchi_id = s.id
  FROM s
 WHERE s.smiles = non.smiles;
commit;


