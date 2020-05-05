/*
 * Create Fragmentation o_source_edge Table SQL Statements :
 *
 * Author | Date    | Version
 * Duncan | 05/2020 | Initial Version
 *
 */

drop table o_source_edge;
create TABLE public.o_source_edge (
    parent_id integer,
    child_id integer,
    label text NOT NULL,
    PRIMARY KEY (parent_id, child_id, label)
);
