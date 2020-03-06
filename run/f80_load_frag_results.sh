#!/bin/bash
# 
# Load Fragmentation results: 
# Purpose: Load Nodes and Edges from Fragmentation step into ISO database 
#
# Parameters:
#    - See file: fragparam.sh and fragpass for fragmentation configuration.
#
# Author | Date    | Version
# Duncan | 03/2020 | Initial Version
#


source fragparam.sh
echo $PYTHONPATH/$FRAGBASEDIR
echo $FRAGNODEFILE
echo $FRAGEDGEFILE

echo "Load table i_node from nodes CSV file"

#\COPY i_node(smiles, hac, rac, child_count, edge_count, ftime) FROM '/data/xchem/nodes.csv' DELIMITER ',' CSV;

psql \
    -X \
    -U postgres \
    -h $DBHOST \
    --echo-all \
    --set AUTOCOMMIT=on \
    --set ON_ERROR_STOP=on \
    -c "\COPY i_node(smiles, hac, rac, child_count, edge_count) FROM '$PYTHONPATH/$FRAGBASEDIR/$FRAGNODEFILE' DELIMITER ',' CSV;" \
    $DATABASE

if [ $? -ne 0 ]; then
    echo "Load Nodes File failed, fault:" 1>&2
    exit $?
fi

echo "Load ISO database from i_nodes"

psql \
    -X \
    -U postgres \
    -h $DBHOST \
    -f f80_load_frag_nodes.sql \
    --echo-all \
    --set AUTOCOMMIT=off \
    --set ON_ERROR_STOP=on \
    $DATABASE

if [ $? -ne 0 ]; then
    echo "Load Nodes results failed, fault:" 1>&2
    exit $?
fi


echo "Load table i_edge from edges CSV file"

#\COPY i_edge(p_smiles, c_smiles, label) FROM '/data/xchem/edges.csv' DELIMITER ',' CSV;

psql \
    -X \
    -U postgres \
    -h $DBHOST \
    --echo-all \
    --set AUTOCOMMIT=on \
    --set ON_ERROR_STOP=on \
    -c "\COPY i_edge(p_smiles, c_smiles, label) FROM '$PYTHONPATH/$FRAGBASEDIR/$FRAGEDGEFILE' DELIMITER ',' CSV;" \
    $DATABASE

if [ $? -ne 0 ]; then
    echo "Load Edges File failed, fault:" 1>&2
    exit $?
fi

echo "Load ISO database from i_edges"

psql \
    -X \
    -U postgres \
    -h $DBHOST \
    -f f80_load_frag_edges.sql \
    --echo-all \
    --set AUTOCOMMIT=off \
    --set ON_ERROR_STOP=on \
    $DATABASE

if [ $? -ne 0 ]; then
    echo "Load Edges results failed, fault:" 1>&2
    exit $?
fi

echo "Load Fragmentation results Successful"
exit 0

