import sys

import duckdb

def count_svs(con, filter_id):
    sql = ''' SELECT svtype, min_svlen, max_svlen
    FROM sv_filters
    WHERE id = ?;
    '''
    filters = con.execute(sql, filter_id).fetchall()

    sql = (f'CREATE TABLE sv_counts_{filter_id}'
           ' FROM'
           ' (SELECT sample, svtype, COUNT(*) AS count'
           ' FROM jrc.joined_raw_calls_svlens'
           ' GROUP BY sample, svtype'
           ' WHERE svtype = ? AND svlen >= ? AND svlen <= ?);'
    con.execute(sql, list(filters[0]))

counts_db = sys.argv[1]
calls_db = sys.argv[2]

with duckdb.connect(counts_db) as con:
    con.execute('ATTACH ? AS jrc;', calls_db)
    filter_ids = con.sql('SELECT DISTINCT(id) FROM sv_filters;').fetchall()
    for i in filter_ids:
        count_svs(con, i[0])
