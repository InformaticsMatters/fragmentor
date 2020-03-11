psql \
    -X \
    -U postgres \
    -h $DBHOST \
    --echo-all \
    --set AUTOCOMMIT=on \
    --set ON_ERROR_STOP=on \
    -c "\COPY i_mols_chemspace(osmiles, isosmiles, nonisosmiles, hac, cmpd_id, price) FROM '$REPPATH/$STANDOUTPUTDIR/$STANDOUTPUTFILE' CSV DELIMITER E'\t' HEADER;" \
    $DATABASE
