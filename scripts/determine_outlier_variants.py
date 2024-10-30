import sys
from pathlib import Path

import duckdb


def find_outliers_from_filter(con, filter_id, min_prop):
    sql = "SELECT svtype FROM sv_filters WHERE id = ?;"
    filters = con.execute(sql, [filter_id]).fetchall()[0]

    sql = (
        f"CREATE OR REPLACE TABLE var_db.outliers_{filter_id}"
        " AS SELECT vid FROM ("
        "SELECT l.vid AS vid,"
        " count(*) AS n_samples,"
        " count(*) FILTER (r.sample IS NOT NULL) AS n_outliers"
        " FROM ("
        "SELECT vid, sample"
        " FROM var_db.variants"
        " WHERE svtype = ?) l"
        f" LEFT JOIN outliers_{filter_id} r"
        " USING(sample)"
        " GROUP BY l.vid)"
        " WHERE n_outliers / n_samples >= ?;"
    )
    con.execute(sql, [filters[0], min_prop])

    sql = (
        f"SELECT DISTINCT l.vid"
        " FROM jrc_db.jrc_clusters l"
        f" JOIN var_db.outliers_{filter_id} r ON (l.member = r.vid)"
    )
    outliers = set([x[0] for x in con.sql(sql).fetchall()])

    return outliers


def find_outlier_variants(con, min_prop):
    sql = "SELECT id FROM sv_filters;"
    filter_ids = con.sql(sql).fetchall()
    count_outliers = [
        find_outliers_from_filter(con, i[0], min_prop) for i in filter_ids
    ]

    sql = (
        "CREATE OR REPLACE TABLE var_db.wgd_outliers"
        " AS SELECT vid FROM ("
        "SELECT l.vid AS vid,"
        " count(*) AS n_samples,"
        " count(*) FILTER (r.sample IS NOT NULL) AS n_outliers"
        " FROM var_db.variants l"
        " LEFT JOIN wgd_outliers r"
        " USING(sample)"
        " GROUP BY l.vid)"
        " WHERE n_outliers / n_samples >= ?;"
    )
    con.execute(sql, [min_prop])

    sql = (
        f"SELECT DISTINCT l.vid"
        " FROM jrc_db.jrc_clusters l"
        f" JOIN var_db.wgd_outliers r ON (l.member = r.vid)"
    )
    wgd_outliers = set([x[0] for x in con.sql(sql).fetchall()])

    return set().union(*count_outliers).union(wgd_outliers)


counts_db = Path(sys.argv[1])
jrc_clusters_db = Path(sys.argv[2])
variants_db_dir = Path(sys.argv[3])
min_outlier_sample_prop = float(sys.argv[4])
if not counts_db.is_file():
    raise FileNotFoundError("Counts database not found")
if not jrc_clusters_db.is_file():
    raise FileNotFoundError("Joined raw calls clusters database not found")
if not variants_db_dir.is_dir():
    raise FileNotFoundError("Variants database directory not found")
if min_outlier_sample_prop < 0 or min_outlier_sample_prop > 1:
    raise ValueError("Min outlier sample proportion must be [0, 1]")


outliers = []
with duckdb.connect(counts_db) as con:
    con.sql(f"ATTACH '{jrc_clusters_db}' AS jrc_db;")
    for db in variants_db_dir.glob("*.duckdb"):
        con.sql(f"ATTACH '{db}' AS var_db;")
        outliers.append(find_outlier_variants(con, min_outlier_sample_prop))
        con.sql("DETACH var_db;")

for vid in set().union(*outliers):
    print(vid)
