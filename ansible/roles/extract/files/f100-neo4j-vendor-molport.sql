/*
 * Extract Neo4j Vendor specific supporting information for source id:
 * Purpose:
 *     Fills the following files with the relevant information from the database
 *     supmolsupedgefile: suppliermol-supplier-edges.csv
 *
 */

--
-- suppliermol-supplier-edges.csv
-- :START_ID(SM_MP),quantity,price_min,price_max,currency,lead_time,:END_ID(S),:TYPE
--
COPY (select code, p.quantity_mg, p.price_min,p.price_max, v.currency, ms.lead_time, v.vendor_name, 'Availability'
        from mol_source ms
        join source s on s.id = ms.source_id
        join vendor_name v on s.name = v.vendor_name
        join price p on molsource_id = ms.id
       where ms.source_id = %(SOURCEID)s)
          to %(SUPMOLSUPEDGEFILE)s DELIMITER ',' CSV;

