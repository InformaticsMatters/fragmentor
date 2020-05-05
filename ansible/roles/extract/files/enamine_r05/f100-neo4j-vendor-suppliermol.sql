/*
 * Extract Neo4j Vendor specific supporting information for source id:
 * Purpose:
 *     Fills the following files with the relevant information from the database
 *     supmolsupedgefile: suppliermol-supplier-edges.csv
 *
 */

--
-- suppliermol-supplier-edges.csv
-- :START_ID(SM),:END_ID(S),:TYPE
--
COPY (select code, v.supplier_node_name, 'Availability'
        from mol_source ms
        join source s on s.id = ms.source_id
        join vendor_name v on s.name = v.vendor_name
       where ms.source_id = %(SOURCEID)s)
          to %(SUPMOLSUPEDGEFILE)s DELIMITER ',' CSV;

