/*
 * Build Neo4j Node Index Data Preparation:
 * Note that this replaces the edge index statements
 * 1. Create parent table o_node_parent if it doesnt exist.
 *    Used as a route for o_node_index_<src> tables to allow combinations
 * 2. Function fills a table with all the parent and child ids from the edge tables
 * that will be contained in the output library for a particular source_id.
 * It is used as a basis for the edge and node-related sql statements.
 *  The edges to be selected will be all the records in the edge table with
 *   a parent in this index.
 *  The nodes to be selected will be all the nodes records with a nonisomol
 *  in this table plus the extra records from mol_source where fragmentation produced no children.
 */

create TABLE IF NOT EXISTS public.o_node_parent (
    nonisomol_id integer,
    PRIMARY KEY (nonisomol_id)
);

CREATE OR REPLACE PROCEDURE extract_o_node_index_src(
   src_id INTEGER,
   chunk_size INTEGER,
   run_offset INTEGER,
   run_limit INTEGER,
   node_table REGCLASS
)
LANGUAGE plpgsql
AS $nodes$
DECLARE
    chunks_count INTEGER := 0;
    commits INTEGER := 0;
    chunks INTEGER := 0;
    src_count INTEGER := 0;
    commit_rate INTEGER := 1;
BEGIN
    chunks := 1 + run_limit / chunk_size;
    RAISE NOTICE 'offset: % has chunks: %', run_offset, chunks;
	for chunk_nr IN 0 .. chunks - 1
    loop
       chunks_count := chunks_count + 1;
       execute 'WITH RECURSIVE source_nodes AS (
           (SELECT nonisomol_id from o_existing_source_mols ms
               limit $2 offset $3)
              union
             select distinct(child_id) as nonisomol_id
             from edge c
             inner join source_nodes p on c.parent_id = p.nonisomol_id
           )
           insert into ' || node_table || '
                select sn.nonisomol_id
                  from source_nodes sn
                order by 1
           ON CONFLICT (nonisomol_id) DO NOTHING;' using src_id, chunk_size, (run_offset+(chunk_size*chunk_nr));
       if chunks_count >= commit_rate then
          chunks_count := 0;
          commit;
          commits := commits + 1;
          RAISE NOTICE 'offset: %, molecules processed: %' , run_offset, commits * commit_rate * chunk_size;
        end if;
    END LOOP;
end $nodes$ ;

