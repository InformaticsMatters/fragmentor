# Fragmentation Enamine - Swich from edges to nodes

## Step 1: Create start point for index extraction

Create o_existing_source_mols table with all nonisomols in mol_source for source_d = 2 that have child_count not null (that is – everything that was successfully fragemented whether they had children or not). 

Enamine tranch 4 is source_id 15.

```
drop table o_existing_source_mols; 
```

```
create table o_existing_source_mols as 
(select distinct m.nonisomol_id from mol_source m 
   join nonisomol n on n.id = m.nonisomol_id and n.child_count notnull 
  where m.source_id = 15);
```

```
ALTER TABLE o_existing_source_mols 
  ADD CONSTRAINT o_existing_source_mols_pk  
    PRIMARY KEY (nonisomol_id);
```

What it does is allow us to: 
a) Speed up the sql statement building each index chunk slightly by removing a subquery (small gain but 1000's times)
b) Add the mol_source child_count = 0 records in directly (simplifying the nodes extract) 

## Step 2: Create table for index extraction

Note because the existing solution merges different vendors, this table is 
implemented as parent table o_edge_parent with child tables o_edge_<vendor>.
Now we are moving to nodes-based index rather than edges this solution will no
longer be used. It needs two fixes:
1. Base the node index on the source_id
2. Combines should now ALL be done in site-combine (see shortcut stories).

So anyway, here I'm just creating the table for the index for source 15:

```
create TABLE public.o_node_index_15 (
    nonisomol_id integer NOT NULL,
    CONSTRAINT "o_node_index_15_pkey" PRIMARY KEY (nonisomol_id)
);
```

## Step 3: Extract node index

Insert/replace procedure in f100-neo4j-source-nodes.sql.
This replaces the edge-related procedure in p10_create_frag_database.sql. 

We have source 39 232 833 nonisomol_id's to gather information for. 
Split over 8 groups of 5 million.
I'm going to use a chunk size of 2000 - this will allow for some big moleculues 
although it will slow the process down a bit.

```
call extract_o_node_index_src(2, 2000, 0, 5000000, 'o_node_index_15' );  
call extract_o_node_index_src(2, 2000, 5000000, 5000000, 'o_node_index_15' );  
call extract_o_node_index_src(2, 2000, 10000000, 5000000, 'o_node_index_15' );  
call extract_o_node_index_src(2, 2000, 15000000, 5000000, 'o_node_index_15' );  
call extract_o_node_index_src(2, 2000, 20000000, 5000000, 'o_node_index_15' );  
call extract_o_node_index_src(2, 2000, 25000000, 5000000, 'o_node_index_15' );  
call extract_o_node_index_src(2, 2000, 30000000, 5000000, 'o_node_index_15' );  
call extract_o_node_index_src(2, 2000, 35000000, 5000000, 'o_node_index_15' );  
```
2501 chunks each.

What it does is allow us to: 
a) Massively increase speed  (nodes are 1 order of magnitude less)
b) Simplify subsequent sql statements  

*Update:* I managed to sucecssfully double the number of jobs and increase the 
chunksize to 3000. These improved the source molecule process rate by about 50%.  
