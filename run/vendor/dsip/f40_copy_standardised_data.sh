psql \
    -X \
    -U postgres \
    -h $DBHOST \
    --echo-all \
    --set AUTOCOMMIT=on \
    --set ON_ERROR_STOP=on \
    -c "\COPY i_mols_dsip(osmiles, isosmiles, nonisosmiles, hac, cmpd_id) FROM '$PYTHONPATH/$STANDOUTPUTDIR/$STANDOUTPUTFILE' CSV DELIMITER E'\t' HEADER;" \
    $DATABASE

