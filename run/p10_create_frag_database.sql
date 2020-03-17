/*
 * Create Fragmentation  Database SQL Statements: 
 * Purpose: Creates a Fragementation  database from scratch and reinitialises information
 *
 * Create Statements originally generated from fairmolecules and dumped 
 * from database version 12.2 by pg_dump version 12.2 (Ubuntu 12.2-2.pgdg18.04+1)
 *
 * Called from p10_create_frag_database.sh
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

DROP TABLE IF EXISTS EDGE;           

DROP TABLE IF EXISTS price;             

DROP TABLE IF EXISTS mol_source;        

DROP TABLE IF EXISTS isomol;            

DROP TABLE IF EXISTS nonisomol;         

DROP TABLE IF EXISTS inchi;             

DROP TABLE IF EXISTS source;            

/*
 * Regenerate database 
 */

--
-- Name: edge; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.edge (
    id integer NOT NULL,
    parent_id integer,
    child_id integer,
    label text NOT NULL
);

--
-- Name: edge_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.edge_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;

--
-- Name: edge_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.edge_id_seq OWNED BY public.edge.id;


--
-- Name: inchi; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.inchi (
    id integer NOT NULL,
    inchik text NOT NULL,
    inchis text NOT NULL
);

--
-- Name: inchi_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.inchi_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;

--
-- Name: inchi_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.inchi_id_seq OWNED BY public.inchi.id;


--
-- Name: isomol; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.isomol (
    id integer NOT NULL,
    smiles text NOT NULL,
    inchik text,
    inchis text,
    nonisomol_id integer NOT NULL
);

--
-- Name: isomol_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.isomol_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;

--
-- Name: isomol_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.isomol_id_seq OWNED BY public.isomol.id;


--
-- Name: mol_source; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.mol_source (
    id integer NOT NULL,
    smiles text NOT NULL,
    code text NOT NULL,
    source_id integer NOT NULL,
    nonisomol_id integer,
    isomol_id integer
);

--
-- Name: mol_source_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.mol_source_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;

--
-- Name: mol_source_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.mol_source_id_seq OWNED BY public.mol_source.id;

--
-- Name: nonisomol; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.nonisomol (
    id integer NOT NULL,
    smiles text NOT NULL,
    inchik text,
    inchis text,
    hac smallint,
    rac smallint,
    child_count smallint,
    edge_count smallint,
    ftime integer,
    inchi_id integer
);

--
-- Name: nonisomol_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.nonisomol_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;

--
-- Name: nonisomol_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.nonisomol_id_seq OWNED BY public.nonisomol.id;


--
-- Name: price; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.price (
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

CREATE SEQUENCE public.price_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;

--
-- Name: price_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.price_id_seq OWNED BY public.price.id;


--
-- Name: source; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.source (
    id integer NOT NULL,
    name text NOT NULL,
    version text NOT NULL,
    currency text
);

--
-- Name: source_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.source_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;

--
-- Name: source_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.source_id_seq OWNED BY public.source.id;

--
-- Name: edge id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.edge ALTER COLUMN id SET DEFAULT nextval('public.edge_id_seq'::regclass);


--
-- Name: inchi id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.inchi ALTER COLUMN id SET DEFAULT nextval('public.inchi_id_seq'::regclass);


--
-- Name: isomol id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.isomol ALTER COLUMN id SET DEFAULT nextval('public.isomol_id_seq'::regclass);


--
-- Name: mol_source id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.mol_source ALTER COLUMN id SET DEFAULT nextval('public.mol_source_id_seq'::regclass);


--
-- Name: nonisomol id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.nonisomol ALTER COLUMN id SET DEFAULT nextval('public.nonisomol_id_seq'::regclass);


--
-- Name: price id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.price ALTER COLUMN id SET DEFAULT nextval('public.price_id_seq'::regclass);


--
-- Name: source id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.source ALTER COLUMN id SET DEFAULT nextval('public.source_id_seq'::regclass);


--
-- Name: edge edge_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.edge
    ADD CONSTRAINT edge_pkey PRIMARY KEY (id);


--
-- Name: inchi inchi_inchis_key; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.inchi
    ADD CONSTRAINT inchi_inchis_key UNIQUE (inchis);


--
-- Name: inchi inchi_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.inchi
    ADD CONSTRAINT inchi_pkey PRIMARY KEY (id);


--
-- Name: isomol isomol_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.isomol
    ADD CONSTRAINT isomol_pkey PRIMARY KEY (id);


--
-- Name: isomol isomol_smiles_key; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.isomol
    ADD CONSTRAINT isomol_smiles_key UNIQUE (smiles);


--
-- Name: mol_source mol_source_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.mol_source
    ADD CONSTRAINT mol_source_pkey PRIMARY KEY (id);


--
-- Name: nonisomol nonisomol_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.nonisomol
    ADD CONSTRAINT nonisomol_pkey PRIMARY KEY (id);


--
-- Name: nonisomol nonisomol_smiles_key; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.nonisomol
    ADD CONSTRAINT nonisomol_smiles_key UNIQUE (smiles);


--
-- Name: price price_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.price
    ADD CONSTRAINT price_pkey PRIMARY KEY (id);


--
-- Name: source source_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.source
    ADD CONSTRAINT source_pkey PRIMARY KEY (id);


--
-- Name: source uq_source; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.source
    ADD CONSTRAINT uq_source UNIQUE (name, version);


--
-- Name: ix_edge_parent_id; Type: INDEX; Schema: public; Owner: postgres
--

CREATE INDEX ix_edge_parent_id ON public.edge USING btree (parent_id);


--
-- Name: ix_inchik; Type: INDEX; Schema: public; Owner: postgres
--

CREATE INDEX ix_inchik ON public.inchi USING hash (inchik);


--
-- Name: ix_isomol_nonisomol_id; Type: INDEX; Schema: public; Owner: postgres
--

CREATE INDEX ix_isomol_nonisomol_id ON public.isomol USING btree (nonisomol_id);


--
-- Name: ix_mol_source_isomol_id; Type: INDEX; Schema: public; Owner: postgres
--

CREATE INDEX ix_mol_source_isomol_id ON public.mol_source USING btree (isomol_id);


--
-- Name: ix_mol_source_nonisomol_id; Type: INDEX; Schema: public; Owner: postgres
--

CREATE INDEX ix_mol_source_nonisomol_id ON public.mol_source USING btree (nonisomol_id);


--
-- Name: ix_mol_source_source_id; Type: INDEX; Schema: public; Owner: postgres
--

CREATE INDEX ix_mol_source_source_id ON public.mol_source USING btree (source_id);


--
-- Name: edge edge_child_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.edge
    ADD CONSTRAINT edge_child_id_fkey FOREIGN KEY (child_id) REFERENCES public.nonisomol(id);


--
-- Name: edge edge_parent_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.edge
    ADD CONSTRAINT edge_parent_id_fkey FOREIGN KEY (parent_id) REFERENCES public.nonisomol(id);


--
-- Name: isomol isomol_nonisomol_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.isomol
    ADD CONSTRAINT isomol_nonisomol_id_fkey FOREIGN KEY (nonisomol_id) REFERENCES public.nonisomol(id);


--
-- Name: mol_source mol_source_isomol_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.mol_source
    ADD CONSTRAINT mol_source_isomol_id_fkey FOREIGN KEY (isomol_id) REFERENCES public.isomol(id);


--
-- Name: mol_source mol_source_nonisomol_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.mol_source
    ADD CONSTRAINT mol_source_nonisomol_id_fkey FOREIGN KEY (nonisomol_id) REFERENCES public.nonisomol(id);


--
-- Name: mol_source mol_source_source_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.mol_source
    ADD CONSTRAINT mol_source_source_id_fkey FOREIGN KEY (source_id) REFERENCES public.source(id) ON DELETE CASCADE;


--
-- Name: nonisomol nonisomol_inchi_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.nonisomol
    ADD CONSTRAINT nonisomol_inchi_id_fkey FOREIGN KEY (inchi_id) REFERENCES public.inchi(id);


--
-- Name: price price_molsource_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.price
    ADD CONSTRAINT price_molsource_id_fkey FOREIGN KEY (molsource_id) REFERENCES public.mol_source(id) ON DELETE CASCADE;

--
-- Table: Source data load.
--
-- Note that these values are fixed as they relate to the vendor parameter numbers.
--

INSERT INTO source (id, name, version, currency)
VALUES (1, 'chembl', '25','');
INSERT INTO source (id, name, version, currency)
VALUES (2, 'chemspace-bb', '2019-07','');
INSERT INTO source (id, name, version, currency)
VALUES (3, 'chemspace-bb', '2019-12','USD');
INSERT INTO source (id, name, version, currency)
VALUES (4, 'molport', '2019-12','USD');
INSERT INTO source (id, name, version, currency)
VALUES (5, 'molport', '2020-01','USD');
INSERT INTO source (id, name, version, currency)
VALUES (6, 'dsip', '1','');

SELECT * FROM source;
