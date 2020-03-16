#!/bin/bash

psql \
    -X \
    -U postgres \
    -h $DBHOST \
    --echo-all \
    --set AUTOCOMMIT=on \
    --set ON_ERROR_STOP=on \
    -c "\COPY i_mols_molport(osmiles, isosmiles, nonisosmiles, hac, cmpd_id, price_1, price_5, price_50, lead_time) FROM '$REPPATH/$STANDOUTPUTDIR/standardised$f.tab' CSV DELIMITER E'\t' HEADER;" \
    $DATABASE