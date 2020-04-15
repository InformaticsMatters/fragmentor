#!/bin/bash
# 
# EXtract Neo4j file:
# Purpose: Extract Node and Edge Information from Fragmentation step into Neo4j files.
#
# Parameters:
#    - See file: fragparam.sh and fragpass for fragmentation configuration.
#
# Author | Date    | Version
# Duncan | 04/2020 | Initial Version
#


source fragparam.sh

echo $FRAGPATH/fragment
echo $NEONODEFILE
echo $NEOEDGEFILE
source $REPPATH/$VENDORPATH/vendorparam.sh

echo "Extracting Neo4j Information Starting"
TSTART=$(date +"%T")
echo "Current time : $TSTART"

echo "Create Neo4jNodes Extract file"

psql \
    -X \
    -U postgres \
    -h $DBHOST \
    --echo-all \
    --set AUTOCOMMIT=off \
    --set ON_ERROR_STOP=on \
    -c "\COPY (WITH RECURSIVE fragments AS ( \
                select parent_id, child_id, parent_smiles, child_smiles, hac, rac, ring_smiles \
                from v_edge \
                where parent_id in ( \
                select n.id from nonisomol n \
                 inner join mol_source ms on ms.nonisomol_id = n.id and ms.source_id = $SOURCEID ) \
                union \
                 select c.parent_id, c.child_id, c.parent_smiles, c.child_smiles, c.hac, c.rac, c.ring_smiles \
                 from v_edge c \
                 inner join fragments p on c.parent_id = p.child_id \
                ) select f1.parent_smiles, f1.hac, f1.rac, f1.ring_smiles, NULL,'F2' as label \
                 from fragments f1 \
                 inner join nonisomol non on non.smiles = f1.parent_smiles \
                union \
                select f2.child_smiles, non.hac, non.rac, non.ring_smiles, NULL, 'F2' \
                 from fragments f2 \
                 inner join nonisomol non on non.smiles = f2.child_smiles \
                 where NOT EXISTS (SELECT 1 FROM v_edge ep WHERE ep.parent_smiles = f2.child_smiles)) \
              TO '$FRAGDATA/fragment/$NEONODEFILE' DELIMITER ',' CSV;"\
    $DATABASE

if [ $? -ne 0 ]; then
    echo "Create Neo4jNodes Extract file failed, fault:" 1>&2
    exit $?
fi

echo "Create Neo4jEdges extract file starting"
TSTART=$(date +"%T")
echo "Current time : $TSTART"


psql \
    -X \
    -U postgres \
    -h $DBHOST \
    --echo-all \
    --set AUTOCOMMIT=off \
    --set ON_ERROR_STOP=on \
    -c "\COPY (with RECURSIVE fragments AS ( \
             select parent_smiles, child_smiles, label, 'FRAG', parent_id, child_id \
             from v_edge \
             where parent_id in ( \
              select n.id from nonisomol n \
               inner join mol_source ms on ms.nonisomol_id = n.id \
               where ms.source_id = $SOURCEID) \
             union \
               select c.parent_smiles, c.child_smiles, c.label, 'FRAG', c.parent_id, c.child_id \
               from v_edge c \
               inner join fragments p on c.parent_id = p.child_id \
            ) select parent_smiles, child_smiles, label, 'FRAG' from fragments) \
              TO '$FRAGDATA/fragment/$NEOEDGEFILE' DELIMITER ',' CSV;" \
    $DATABASE

if [ $? -ne 0 ]; then
    echo "Create Neo4jEdges extract file failed, fault:" 1>&2
    exit $?
fi

echo "Extracting Neo4j Information Successful"
TEND=$(date +"%T")
echo "Current time : $TEND"

exit 0

