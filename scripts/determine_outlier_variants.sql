.mode tabs
.headers off
SELECT vid
FROM (
    SELECT l.vid AS vid, COUNT(*) FILTER(r.sample IS NULL) AS nils
    FROM variants l
    LEFT JOIN outlier_samples r
    USING (sample)
    GROUP BY vid
) 
WHERE nils = 0;
