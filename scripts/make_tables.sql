CREATE TABLE joined_raw_calls_svlens (
    vid VARCHAR,
    svtype VARCHAR,
    svlen INTEGER,
    sample VARCHAR
);

CREATE TABLE joined_raw_calls_clusters (
    vid VARCHAR,
    member VARCHAR
);

COPY joined_raw_calls_svlens
FROM 'joined_raw_calls_svlens.tsv' (
    FORMAT CSV,
    DELIMITER '\t',
    HEADER false
);

COPY joined_raw_calls_clusters
FROM 'joined_raw_calls_clusters.tsv' (
    FORMAT CSV,
    DELIMITER '\t',
    HEADER false
);
