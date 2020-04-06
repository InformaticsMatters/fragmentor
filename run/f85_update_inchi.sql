/*
 * Load Fragmented Nodes Results SQL Statements: 
 * Purpose: Loads fragmented node data into Frag database
 *
 * Called from p80_load_nodes.sh
 *
 * Author | Date    | Version
 * Duncan | 03/2020 | Initial Version
 *
 */

\timing

/*
 * Insert INCHI tables.
 */
begin;
SELECT clock_timestamp();

INSERT INTO inchi (inchik, inchis)
  SELECT sinchik, sinchis from i_inchi
  ON CONFLICT ON CONSTRAINT inchi_inchis_key DO NOTHING;
commit;

/*
 * Update inchi key for nonisomols inserted in f40.
 */
begin;
SELECT clock_timestamp();

WITH s
    AS (SELECT i.smiles, i.ninchik, i.ninchis, inc.id
  FROM i_inchi i
  JOIN inchi inc ON inc.inchik = i.sinchik
 WHERE i.ninchik <> '')
UPDATE nonisomol non
   SET inchik = s.ninchik, inchis = s.ninchis, inchi_id = s.id
  FROM s
 WHERE s.smiles = non.smiles;
commit;


/*
 * Update isomol from i_node with the inchi key.
 */
begin;
SELECT clock_timestamp();

WITH s
    AS (SELECT i.smiles, i.ninchik, i.ninchis
  FROM i_inchi i
  JOIN isomol iso ON iso.smiles = i.smiles
 WHERE i.ninchik <> '')
UPDATE isomol iso
   SET inchik = s.ninchik, inchis = s.ninchis
  FROM s
 WHERE s.smiles = iso.smiles;
commit;
