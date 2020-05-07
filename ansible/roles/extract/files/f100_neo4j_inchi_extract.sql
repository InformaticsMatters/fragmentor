/*
 * Extract Neo4j Inchi Data SQL Statements:
 * Purpose: Extract Neo4j Inchi into CSV file on database server
 */

--
-- inchi.csv
-- inchid:ID(Inchi),inchik,inchis,:LABEL
-- Inchis used by molecules for source + label "Inchi"
---

COPY (select np.inchi_id, inc1.inchik, inc1.inchis, 'Inchi' as label
       from o_source_edge os
       join nonisomol np ON np.id = os.parent_id
       join inchi inc1 on inc1.id = np.inchi_id
       left join mol_source ms on os.parent_id = ms.nonisomol_id and ms.source_id = %(SOURCEID)s
      union
     select nc.inchi_id, inc2.inchik, inc2.inchis, 'Inchi'
       from o_source_edge os2
       join nonisomol nc ON nc.id = os2.child_id and nc.child_count = 0
       join inchi inc2 on inc2.id = nc.inchi_id)
   TO %(INCHIIDFILE)s DELIMITER ',' CSV;
