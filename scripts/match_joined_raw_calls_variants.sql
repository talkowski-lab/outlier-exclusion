CREATE TEMP TABLE t1 (
    vid VARCHAR
);

COPY t1 
FROM 'outlier_variants.list' (
    FORMAT CSV,
    DELIMITER '\t',
    HEADER false
);

.mode tabs
SELECT jrc.vid AS joined_raw_calls_vid
FROM joined_raw_calls_svlens jrc
JOIN t1 ON (jrc.member = t1.vid);
