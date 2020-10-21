/*
 * Load configuration data SQL Statements:
 * Purpose: Initialises venfor information
 *
 * Author | Date    | Version
 * Duncan | 05/2020 | Initial Version
 *
 */

insert into vendor_name (vendor_name, currency, supplier_node_name, supplier_node_label) values ('xchem_dsip', NULL, 'Xchem-Dsip','V_XDSIP');
insert into vendor_name (vendor_name, currency, supplier_node_name, supplier_node_label) values ('chemspace_bb', 'USD', 'ChemSpace-BB','V_CS_BB');
insert into vendor_name (vendor_name, currency, supplier_node_name, supplier_node_label) values ('molport', 'USD', 'MolPort','V_MP');
insert into vendor_name (vendor_name, currency, supplier_node_name, supplier_node_label) values ('enamine_ro5', NULL, 'REAL','V_REAL');
insert into vendor_name (vendor_name, currency, supplier_node_name, supplier_node_label) values ('xchem_spot', NULL, 'Xchem-Spot','V_XSPOT');
insert into vendor_name (vendor_name, currency, supplier_node_name, supplier_node_label) values ('xchem_probe', NULL, 'Xchem-Probe','V_XPROB');
