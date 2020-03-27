/*
 * Load Fragmented Results SQL Statements: 
 * Purpose: Loads fragmented node and edge data into Frag database
 *
 * Called from p80_load_frag_results.sh
 *
 * Author | Date    | Version
 * Duncan | 03/2020 | Initial Version
 *
 */

/*
 * Check which edges are already present in the database 
 */

\timing

begin;

/*
 * Create i_edge_present
 */
SELECT clock_timestamp();

DROP TABLE IF EXISTS i_edge;
CREATE TABLE i_edge_present (
  p_smiles TEXT,
  c_smiles TEXT,
  label TEXT
);

/*
 * Delete rows from i_edge that already exist in edge and insert them into i_edge_present.
 */

WITH presnt AS (DELETE FROM i_edge i RETURNING p_smiles, c_smiles, label
              WHERE EXISTS
                (SELECT 1 FROM edge e
                   JOIN nonisomol np ON np.id = e.parent_id
                   JOIN nonisomol nc ON nc.id = e.child_id
                  WHERE i.p_smiles = np.smiles AND i.c_smiles = nc.smiles)
    INSERT INTO i_edge_present (p_smiles, c_smiles, label) SELECT * FROM presnt;
commit;

/*
 * count present edges  
 */
SELECT count(*) FROM i_edge_present;
SELECT count(*) FROM i_edge;
/*
 * Create i_edge index
 */
begin;
SELECT clock_timestamp();

DROP INDEX IF EXISTS ix_i_edge_psmiles;
DROP INDEX IF EXISTS ix_i_edge_csmiles;
CREATE INDEX ix_i_edge_psmiles ON i_edge USING btree (p_smiles);
CREATE INDEX ix_i_edge_csmiles ON i_edge USING btree (c_smiles);
commit;

/*
 * Load new edges 
 */
begin;
SELECT clock_timestamp();

INSERT INTO edge (parent_id, child_id, label)
  SELECT np.id, nc.id, i.label FROM i_edge i
    JOIN nonisomol np ON np.smiles = i.p_smiles 
    JOIN nonisomol nc ON nc.smiles = i.c_smiles;
commit;

