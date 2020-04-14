psql \
    -X \
    -U postgres \
    -h $DBHOST \
    --echo-all \
    --set AUTOCOMMIT=on \
    --set ON_ERROR_STOP=on \
    -c "\COPY i_mols_enamine(osmiles, isosmiles, nonisosmiles, hac, cmpd_id) FROM '$STANDPATH/standardise/standardised$f.tab' CSV DELIMITER E'\t' HEADER;" \
    $DATABASE

