/*
 * Extract Neo4j Supporting information for source id:
 * Purpose:
 *     Fills the following files with the relevant information from the database
 *     supnodefile: supplier-nodes.csv
 *     supmolnodefile: suppliermol-nodes.csv
 *     supmolsupedgefile: suppliermol-supplier-edges.csv
 *     molsupmoledgefile: molecule-suppliermol-edges.csv
 *     isosupmoledgefile: isomol-suppliermol-edges.csv
 *     isomoledgefile: isomol-molecule-edges.csv
 *     isomolnodefile: isomol-nodes.csv
 *
 */

--
-- Supplier Nodes.csv
-- name:ID(S),graph_version,processing_version,process_id,build_number:int,limit:int,min_hac:int,max_hac:int,build_datetime:datetime,label,:LABEL
-- NB - suggest these values are ultimately added to "Vendor_name" and Source.
--
COPY (select v.supplier_node_name, %(GRAPHVERSION)s ,s.name || '/' || s.version, %(PROCESSID)s, 1,
             0, s.min_hac, s.max_hac, s.start_datetime, v.supplier_node_label, 'Supplier'
        from source s
        join vendor_name v on s.name = v.vendor_name
       where s.id = %(SOURCEID)s)
          to %(SUPNODEFILE)s DELIMITER ',' CSV;

--
-- suppliermol-nodes.csv
-- cmpd_id:ID(SM),osmiles,:LABEL
--
COPY (select code, smiles, 'Available'
        from mol_source ms
       where ms.source_id = %(SOURCEID)s)
          to %(SUPMOLNODEFILE)s DELIMITER ',' CSV;

--
-- suppliermol-supplier-edges.csv
-- :START_ID(SM),:END_ID(S),:TYPE
--
COPY (select code, 'Xchem', 'Availability'
        from mol_source ms
       where ms.source_id = %(SOURCEID)s)
          to %(SUPMOLSUPEDGEFILE)s DELIMITER ',' CSV;

--
-- molecule-suppliermol-edges.csv
-- :START_ID(F2),:END_ID(SM),:TYPE
--
COPY (select non.smiles, ms.code, 'HasVendor'
        from mol_source ms
        join nonisomol non on ms.nonisomol_id = non.id
       where ms.isomol_id is null
         and ms.source_id = %(SOURCEID)s)
          to %(MOLSUPMOLEDGEFILE)s DELIMITER ',' CSV;

--
-- isomol-suppliermol-edges.csv
-- :START_ID(ISOMOL),:END_ID(SM),:TYPE
--
COPY (select iso.smiles, ms.code, 'HasVendor'
        from mol_source ms
        join isomol iso on ms.isomol_id = iso.id
         and ms.source_id = %(SOURCEID)s)
          to %(ISOSUPMOLEDGEFILE)s DELIMITER ',' CSV;

--
-- isomol-molecule-edges.csv
-- :START_ID(ISOMOL),:END_ID(F2),:TYPE
--
COPY (select iso.smiles, non.smiles, 'NonIso'
        from mol_source ms
        join isomol iso on iso.id = ms.isomol_id
        join nonisomol non on non.id = iso.nonisomol_id
         and ms.source_id = %(SOURCEID)s)
          to %(ISOMOLEDGEFILE)s DELIMITER ',' CSV;

--
-- isomolnodefile: isomol-nodes.csv
-- smiles:ID(ISOMOL),cmpd_ids:STRING[],:LABEL
--
COPY (select iso.smiles, ms.code, 'CanSmi;Mol;XChem'
        from mol_source ms
        join isomol iso on iso.id = ms.isomol_id
         and ms.source_id = %(SOURCEID)s)
          to %(ISOMOLNODEFILE)s DELIMITER ',' CSV;
