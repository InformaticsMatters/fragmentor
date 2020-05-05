/*
 * Create Fragmentation  Database SQL Statements: 
 * Purpose: Creates a Fragementation  database from scratch and reinitialises information
 *
 * Create Statements originally generated from fairmolecules and dumped 
 * from database version 12.2 by pg_dump version 12.2 (Ubuntu 12.2-2.pgdg18.04+1)
 *
 * NOTE: All tables are created under the public owner. If access permissions 
 * are required, then these will have to be set.
 * 
 * Author | Date    | Version  
 * Duncan | 03/2020 | Initial Version
 *
 */

/*
 * Drop all tables if they exist (NB Order because of table key constrants) 
 */


drop view IF EXISTS v_edge;

drop view IF EXISTS v_edge_node;

drop table IF EXISTS edge;

drop table IF EXISTS price;

drop table IF EXISTS mol_source;

drop table IF EXISTS isomol;

drop table IF EXISTS nonisomol;

drop table IF EXISTS inchi;

drop table IF EXISTS source;

/*
 * Regenerate database 
 */

--
-- Name: edge; Type: TABLE; Schema: public; Owner: postgres
--

create TABLE public.edge (
    id integer NOT NULL,
    parent_id integer,
    child_id integer,
    label text NOT NULL,
    source_id integer
);

--
-- Name: edge_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

create sequence public.edge_id_seq
    AS integer
    start with 1
    increment by 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;

--
-- Name: edge_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

alter sequence public.edge_id_seq OWNED BY public.edge.id;


--
-- Name: inchi; Type: TABLE; Schema: public; Owner: postgres
--

create TABLE public.inchi (
    id integer NOT NULL,
    inchik text NOT NULL,
    inchis text NOT NULL
);

--
-- Name: inchi_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

create sequence public.inchi_id_seq
    AS integer
    start with 1
    increment by 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;

--
-- Name: inchi_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

alter sequence public.inchi_id_seq OWNED BY public.inchi.id;


--
-- Name: isomol; Type: TABLE; Schema: public; Owner: postgres
--

create TABLE public.isomol (
    id integer NOT NULL,
    smiles text NOT NULL,
    inchik text,
    inchis text,
    nonisomol_id integer NOT NULL,
    source_id integer
);

--
-- Name: isomol_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

create sequence public.isomol_id_seq
    AS integer
    start with 1
    increment by 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;

--
-- Name: isomol_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

alter sequence public.isomol_id_seq OWNED BY public.isomol.id;


--
-- Name: mol_source; Type: TABLE; Schema: public; Owner: postgres
--

create TABLE public.mol_source (
    id integer NOT NULL,
    smiles text NOT NULL,
    code text NOT NULL,
    source_id integer NOT NULL,
    lead_time INTEGER,
    nonisomol_id integer,
    isomol_id integer
);

--
-- Name: mol_source_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

create sequence public.mol_source_id_seq
    AS integer
    start with 1
    increment by 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;

--
-- Name: mol_source_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

alter sequence public.mol_source_id_seq OWNED BY public.mol_source.id;

--
-- Name: nonisomol; Type: TABLE; Schema: public; Owner: postgres
--

create TABLE public.nonisomol (
    id integer NOT NULL,
    smiles text NOT NULL,
    inchik text,
    inchis text,
    hac smallint,
    rac smallint,
    ring_smiles text,
    child_count smallint,
    edge_count smallint,
    ftime integer,
    inchi_id integer,
    source_id integer
);

--
-- Name: nonisomol_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

create sequence public.nonisomol_id_seq
    AS integer
    start with 1
    increment by 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;

--
-- Name: nonisomol_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

alter sequence public.nonisomol_id_seq OWNED BY public.nonisomol.id;


--
-- Name: price; Type: TABLE; Schema: public; Owner: postgres
--

create TABLE public.price (
    id integer NOT NULL,
    quantity_mg integer,
    price integer,
    price_min integer,
    price_max integer,
    molsource_id integer
);

--
-- Name: price_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

create sequence public.price_id_seq
    AS integer
    start with 1
    increment by 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;

--
-- Name: price_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

alter sequence public.price_id_seq OWNED BY public.price.id;


--
-- Name: source; Type: TABLE; Schema: public; Owner: postgres
--

create TABLE public.source (
    id integer NOT NULL,
    name text NOT NULL,
    version text NOT NULL,
    frag_limit integer,
    min_hac integer,
    max_hac integer,
    max_frags integer,
    start_datetime timestamp,
    run_summary_log BYTEA
);

--
-- Name: source_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

create sequence public.source_id_seq
    AS integer
    start with 1
    increment by 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;

--
-- Name: source_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

alter sequence public.source_id_seq OWNED BY public.source.id;

--
-- Name: vendor_version; Type: TABLE; Schema: public; Owner: postgres
--

create TABLE public.vendor_name (
    vendor_name text NOT NULL,
    currency text,
    supplier_node_name text,
    supplier_node_label text,
    PRIMARY KEY (vendor_name)
);

--
-- Name: edge id; Type: DEFAULT; Schema: public; Owner: postgres
--

alter table ONLY public.edge alter COLUMN id SET DEFAULT nextval('public.edge_id_seq'::regclass);


--
-- Name: inchi id; Type: DEFAULT; Schema: public; Owner: postgres
--

alter table ONLY public.inchi alter COLUMN id SET DEFAULT nextval('public.inchi_id_seq'::regclass);


--
-- Name: isomol id; Type: DEFAULT; Schema: public; Owner: postgres
--

alter table ONLY public.isomol alter COLUMN id SET DEFAULT nextval('public.isomol_id_seq'::regclass);


--
-- Name: mol_source id; Type: DEFAULT; Schema: public; Owner: postgres
--

alter table ONLY public.mol_source alter COLUMN id SET DEFAULT nextval('public.mol_source_id_seq'::regclass);


--
-- Name: nonisomol id; Type: DEFAULT; Schema: public; Owner: postgres
--

alter table ONLY public.nonisomol alter COLUMN id SET DEFAULT nextval('public.nonisomol_id_seq'::regclass);


--
-- Name: price id; Type: DEFAULT; Schema: public; Owner: postgres
--

alter table ONLY public.price alter COLUMN id SET DEFAULT nextval('public.price_id_seq'::regclass);


--
-- Name: source id; Type: DEFAULT; Schema: public; Owner: postgres
--

alter table ONLY public.source alter COLUMN id SET DEFAULT nextval('public.source_id_seq'::regclass);


--
-- Name: edge edge_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

alter table ONLY public.edge
    ADD CONSTRAINT edge_pkey PRIMARY KEY (id);


--
-- Name: inchi inchi_inchis_key; Type: CONSTRAINT; Schema: public; Owner: postgres
--

alter table ONLY public.inchi
    ADD CONSTRAINT inchi_inchis_key UNIQUE (inchis);


--
-- Name: inchi inchi_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

alter table ONLY public.inchi
    ADD CONSTRAINT inchi_pkey PRIMARY KEY (id);


--
-- Name: isomol isomol_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

alter table ONLY public.isomol
    ADD CONSTRAINT isomol_pkey PRIMARY KEY (id);


--
-- Name: isomol isomol_smiles_key; Type: CONSTRAINT; Schema: public; Owner: postgres
--

alter table ONLY public.isomol
    ADD CONSTRAINT isomol_smiles_key UNIQUE (smiles);


--
-- Name: mol_source mol_source_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

alter table ONLY public.mol_source
    ADD CONSTRAINT mol_source_pkey PRIMARY KEY (id);


--
-- Name: nonisomol nonisomol_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

alter table ONLY public.nonisomol
    ADD CONSTRAINT nonisomol_pkey PRIMARY KEY (id);


--
-- Name: nonisomol nonisomol_smiles_key; Type: CONSTRAINT; Schema: public; Owner: postgres
--

alter table ONLY public.nonisomol
    ADD CONSTRAINT nonisomol_smiles_key UNIQUE (smiles);


--
-- Name: price price_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

alter table ONLY public.price
    ADD CONSTRAINT price_pkey PRIMARY KEY (id);


--
-- Name: source source_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

alter table ONLY public.source
    ADD CONSTRAINT source_pkey PRIMARY KEY (id);


--
-- Name: source uq_source; Type: CONSTRAINT; Schema: public; Owner: postgres
--

alter table ONLY public.source
    ADD CONSTRAINT uq_source UNIQUE (name, version);


--
-- Name: ix_edge_parent_id; Type: INDEX; Schema: public; Owner: postgres
--

create INDEX ix_edge_parent_id ON public.edge USING btree (parent_id);

--
-- Name: ix_edge_parent_child; Type: INDEX; Schema: public; Owner: postgres
--
create INDEX ix_edge_parent_child ON public.edge USING btree (parent_id, child_id);

--
-- Name: ix_inchik; Type: INDEX; Schema: public; Owner: postgres
--

create INDEX ix_inchik ON public.inchi USING hash (inchik);


--
-- Name: ix_isomol_nonisomol_id; Type: INDEX; Schema: public; Owner: postgres
--

create INDEX ix_isomol_nonisomol_id ON public.isomol USING btree (nonisomol_id);


--
-- Name: ix_mol_source_isomol_id; Type: INDEX; Schema: public; Owner: postgres
--

create INDEX ix_mol_source_isomol_id ON public.mol_source USING btree (isomol_id);


--
-- Name: ix_mol_source_nonisomol_id; Type: INDEX; Schema: public; Owner: postgres
--

create INDEX ix_mol_source_nonisomol_id ON public.mol_source USING btree (nonisomol_id);


--
-- Name: ix_mol_source_source_id; Type: INDEX; Schema: public; Owner: postgres
--

create INDEX ix_mol_source_source_id ON public.mol_source USING btree (source_id);


--
-- Name: edge edge_child_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

alter table ONLY public.edge
    ADD CONSTRAINT edge_child_id_fkey FOREIGN KEY (child_id) REFERENCES public.nonisomol(id);


--
-- Name: edge edge_parent_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

alter table ONLY public.edge
    ADD CONSTRAINT edge_parent_id_fkey FOREIGN KEY (parent_id) REFERENCES public.nonisomol(id);


--
-- Name: isomol isomol_nonisomol_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

alter table ONLY public.isomol
    ADD CONSTRAINT isomol_nonisomol_id_fkey FOREIGN KEY (nonisomol_id) REFERENCES public.nonisomol(id);


--
-- Name: mol_source mol_source_isomol_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

alter table ONLY public.mol_source
    ADD CONSTRAINT mol_source_isomol_id_fkey FOREIGN KEY (isomol_id) REFERENCES public.isomol(id);


--
-- Name: mol_source mol_source_nonisomol_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

alter table ONLY public.mol_source
    ADD CONSTRAINT mol_source_nonisomol_id_fkey FOREIGN KEY (nonisomol_id) REFERENCES public.nonisomol(id);


--
-- Name: mol_source mol_source_source_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

alter table ONLY public.mol_source
    ADD CONSTRAINT mol_source_source_id_fkey FOREIGN KEY (source_id) REFERENCES public.source(id) ON delete CASCADE;


--
-- Name: nonisomol nonisomol_inchi_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

alter table ONLY public.nonisomol
    ADD CONSTRAINT nonisomol_inchi_id_fkey FOREIGN KEY (inchi_id) REFERENCES public.inchi(id);


--
-- Name: price price_molsource_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

alter table ONLY public.price
    ADD CONSTRAINT price_molsource_id_fkey FOREIGN KEY (molsource_id) REFERENCES public.mol_source(id) ON delete CASCADE;

--
-- View Used in extracting data for Neo4j database - only used in old recursive extract process
--

create view v_edge as
    SELECT e.id, e.parent_id, e.child_id,
           p.smiles as parent_smiles,
           c.smiles as child_smiles,
           e.label
      FROM edge e
INNER JOIN nonisomol p ON e.parent_id = p.id
INNER JOIN nonisomol c ON e.child_id = c.id;

create view v_edge_node as
    SELECT e.id, e.parent_id, e.child_id,
           p.smiles as parent_smiles,
           c.smiles as child_smiles,
           p.hac as hac,
           p.rac as rac,
           p.ring_smiles as ring_smiles,
           c.hac as child_hac,
           c.rac as child_rac,
           c.ring_smiles as child_ring_smiles,
           c.child_count as child_child_count,
           e.label
      FROM edge e
INNER JOIN nonisomol p ON e.parent_id = p.id
INNER JOIN nonisomol c ON e.child_id = c.id;

--
-- Tuning Parameters
--

-- This is set for the largest tables so that autovacuuming starts more quickly.

alter table nonisomol set (autovacuum_vacuum_scale_factor = 0.01);
alter table nonisomol set (autovacuum_vacuum_cost_limit = 2000);
alter table edge set (autovacuum_vacuum_scale_factor = 0.01);
alter table edge set (autovacuum_vacuum_cost_limit = 2000);
alter table mol_source set (autovacuum_vacuum_scale_factor = 0.01);
alter table mol_source set (autovacuum_vacuum_cost_limit = 2000);


--
-- Stored Procedure for loading o_source_edge
--

CREATE OR REPLACE PROCEDURE extract_o_source_edge(
   src_id INTEGER,
   commit_rate INTEGER,
   chunk_size INTEGER
)
LANGUAGE plpgsql
AS $edges$
DECLARE
    chunks_count INTEGER := 0;
    commits INTEGER := 0;
    chunks INTEGER := 0;
    src_count INTEGER := 0;
BEGIN
	select count(*) into src_count from mol_source where source_id = src_id;
    chunks := 1 + src_count / chunk_size;
    RAISE NOTICE 'chunks %', chunks;
	for chunk_nr IN 0 .. chunks - 1
    loop
       chunks_count := chunks_count + 1;
       WITH RECURSIVE source_edges AS (
           select parent_id, child_id, label
             from edge e
            where e.parent_id in (SELECT ms.nonisomol_id
                                   from mol_source ms
                                   where ms.source_id = src_id limit chunk_size offset (chunk_size*chunk_nr))
            union
           select c.parent_id, c.child_id, c.label
             from edge c
             inner join source_edges p on c.parent_id = p.child_id)
           insert into o_source_edge
              select parent_id, child_id, label
                from source_edges
           ON CONFLICT (parent_id, child_id, label) DO NOTHING;
       if chunks_count >= commit_rate then
          chunks_count := 0;
          commit;
          commits := commits + 1;
          RAISE NOTICE 'molecules processed: %', commits * commit_rate * chunk_size;
       end if;
    END LOOP;
end $edges$ ;




