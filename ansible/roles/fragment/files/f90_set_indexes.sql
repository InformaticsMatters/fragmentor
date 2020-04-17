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

/*
 * Create i_edge index
 */
begin;
SELECT clock_timestamp();
CREATE INDEX ix_i_edge_psmiles ON i_edge USING btree (p_smiles);
CREATE INDEX ix_i_edge_csmiles ON i_edge USING btree (c_smiles);
commit;

