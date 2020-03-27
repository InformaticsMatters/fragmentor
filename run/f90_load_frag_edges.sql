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
/*
 * Create i_edge index
 */
begin;
SELECT clock_timestamp();
CREATE INDEX ix_i_edge_psmiles ON i_edge USING btree (p_smiles);
CREATE INDEX ix_i_edge_csmiles ON i_edge USING btree (c_smiles);
commit;

begin;
SELECT clock_timestamp();

UPDATE i_edge i SET present = TRUE WHERE EXISTS
  (SELECT 1 FROM edge e 
  JOIN nonisomol np ON np.id = e.parent_id 
  JOIN nonisomol nc ON nc.id = e.child_id 
  WHERE i.p_smiles = np.smiles AND i.c_smiles = nc.smiles);

/*
 * count present edges
 */
SELECT present, count(*) FROM i_edge GROUP BY present;

commit;

/*
 * Load new edges 
 */
begin;
SELECT clock_timestamp();

INSERT INTO edge (parent_id, child_id, label)
  SELECT np.id, nc.id, i.label FROM i_edge i
    JOIN nonisomol np ON np.smiles = i.p_smiles 
    JOIN nonisomol nc ON nc.smiles = i.c_smiles
    WHERE i.present IS NULL;
commit;

