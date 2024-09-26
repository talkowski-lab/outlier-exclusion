CREATE TABLE join_raw_calls_svlens (
    vid VARCHAR,
    svtype VARCHAR,
    svlen INTEGER,
    sample VARCHAR
);

CREATE TABLE join_raw_calls_clusters (
    vid VARCHAR,
    member VARCHAR
);

COPY join_raw_calls_svlens
FROM 'join_raw_calls_svlens.tsv' (
    FORMAT CSV,
    DELIMITER '\t',
    HEADER false
);

COPY join_raw_calls_clusters
FROM 'join_raw_calls_clusters.tsv' (
    FORMAT CSV,
    DELIMITER '\t',
    HEADER false
);