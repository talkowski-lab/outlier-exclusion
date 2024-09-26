COPY
(SELECT sample, svtype, COUNT(*) AS count
FROM join_raw_calls_svlens
GROUP BY sample, svtype)
TO 'joined_raw_calls_svcounts_ALL.tsv'
(DELIMITER '\t', HEADER true);

PREPARE query_svtype AS
SELECT sample, COUNT(*) AS count
FROM join_raw_calls_svlens
WHERE svtype = ? AND svlen >= ? AND svlen <= ?
GROUP BY sample;

.mode tabs
.once 'joined_raw_calls_svcount_{{svtype}}.tsv'
EXECUTE query_svtype('{{svtype}}', {{min_svlen}}, {{max_svlen}});
