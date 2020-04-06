#!/bin/bash
# 
# Load Fragmentation results: 
# Purpose: Load Edges from Fragmentation step into Frag database
# 1. Split edges file into chunks of chunk_size EDGECHUNK for loading into the frag database
# 2. For each chunk copy the data into i_edge
# 3. Create indexes on i_edge (this may not be necessary)
# 4. Extract the edge infomation into a newedgechunk* csv file
#
# 5. Drop indexes on edge table for speed of loading
# 6. For each newedgechunkfile, load the edge file.
# 7. Recreate the indexes on the edge table for the next run.
#
# Parameters:
#    - See file: fragparam.sh and fragpass for edgechunk configuration.
#
# Author | Date    | Version
# Duncan | 03/2020 | Initial Version
#


source fragparam.sh
echo $REPPATH/$FRAGBASEDIR
echo $FRAGEDGEFILE

echo "Loading Edges Starting"
TSTART=$(date +"%T")
echo "Current time : $TSTART"

read lines filename <<< $(wc -l $REPPATH/$FRAGBASEDIR/$FRAGEDGEFILE)
echo "lines=$lines filename=$filename"

cat $REPPATH/$FRAGBASEDIR/$FRAGEDGEFILE | split -d -l $EDGECHUNK - $REPPATH/$FRAGBASEDIR/edgechunk_

for f in $REPPATH/$FRAGBASEDIR/edgechunk*; do

    echo "Processing Filename $f"
    TSTART=$(date +"%T")
    echo "Current time : $TSTART"

    echo "Create table i_edge"

    psql \
        -X \
        -U postgres \
        -h $DBHOST \
        -f f90_create_i_edges.sql \
        --echo-all \
        --set AUTOCOMMIT=on \
        --set ON_ERROR_STOP=on \
        $DATABASE

    if [ $? -ne 0 ]; then
        echo "Create table i_edge failed, fault:" 1>&2
        exit $?
    fi

    echo "Load table i_edge from edges CSV file"
    TSTART=$(date +"%T")
    echo "Current time : $TSTART"

    psql \
        -X \
        -U postgres \
        -h $DBHOST \
        --echo-all \
        --set AUTOCOMMIT=on \
        --set ON_ERROR_STOP=on \
        -c "\COPY i_edge(p_smiles, c_smiles, label) FROM '$f' DELIMITER ',' CSV;" \
        $DATABASE

    if [ $? -ne 0 ]; then
        echo "Load Edges File failed, fault:" 1>&2
        exit $?
    fi

     echo "Create indexes"
     TSTART=$(date +"%T")
     echo "Current time : $TSTART"

     psql \
         -X \
         -U postgres \
         -h $DBHOST \
         -f f90_set_indexes.sql \
         --echo-all \
         --set AUTOCOMMIT=off \
         --set ON_ERROR_STOP=on \
         $DATABASE

     if [ $? -ne 0 ]; then
         echo "Create indexes failed, fault:" 1>&2
         exit $?
     fi

     echo "Extract file newedgechunk*.csv From i_edge"
     TSTART=$(date +"%T")
     echo "Current time : $TSTART"

     psql \
        -X \
        -U postgres \
        -h $DBHOST \
        --echo-all \
        --set AUTOCOMMIT=off \
        --set ON_ERROR_STOP=on \
        -c "\COPY (SELECT np.id, nc.id, ied.label \
                    FROM i_edge ied \
                    JOIN nonisomol np ON np.smiles = ied.p_smiles \
                    JOIN nonisomol nc ON nc.smiles = ied.c_smiles \
                   WHERE NOT EXISTS \
                    (SELECT 1 FROM edge edg \
                      WHERE edg.parent_id = np.id  \
                        AND edg.child_id  = nc.id)) \
              TO '$f.new' DELIMITER ',' CSV;"\
        $DATABASE

     if [ $? -ne 0 ]; then
        echo "Extract Edges File failed, fault:" 1>&2
        exit $?
     fi

done

echo "Drop indexes from Edge table"
TSTART=$(date +"%T")
echo "Current time : $TSTART"

psql \
  -X \
  -U postgres \
  -h $DBHOST \
  --echo-all \
  --set AUTOCOMMIT=on \
  --set ON_ERROR_STOP=on \
  -c "DROP INDEX IF EXISTS ix_edge_parent_id;
      DROP INDEX IF EXISTS ix_edge_parent_child;" \
  $DATABASE

if [ $? -ne 0 ]; then
  echo "Drop indexes failed, fault:" 1>&2
  exit $?
fi

echo "Load table edge from New Edges CSV files"


for f in $REPPATH/$FRAGBASEDIR/edgechunk_*.new; do

    echo "Processing Filename $f"
    TSTART=$(date +"%T")
    echo "Current time : $TSTART"

    psql \
        -X \
        -U postgres \
        -h $DBHOST \
        --echo-all \
        --set AUTOCOMMIT=on \
        --set ON_ERROR_STOP=on \
        -c "\COPY edge(parent_id, child_id, label) FROM '$f' DELIMITER ',' CSV;" \
        $DATABASE

    if [ $? -ne 0 ]; then
        echo "Load Edges File failed, fault:" 1>&2
        exit $?
    fi

done

echo "Drop indexes from Edge table"
TSTART=$(date +"%T")
echo "Current time : $TSTART"

psql \
  -X \
  -U postgres \
  -h $DBHOST \
  --echo-all \
  --set AUTOCOMMIT=on \
  --set ON_ERROR_STOP=on \
  -c "CREATE INDEX ix_edge_parent_id ON public.edge USING btree (parent_id);
      CREATE INDEX ix_edge_parent_child ON public.edge USING btree (parent_id, child_id);" \
  $DATABASE

if [ $? -ne 0 ]; then
  echo "Drop indexes failed, fault:" 1>&2
  exit $?
fi

rm $REPPATH/$FRAGBASEDIR/edgechunk_*

echo "Loading Edges Successful"
TEND=$(date +"%T")
echo "Current time : $TEND"

exit 0
