"""Determine GATK-SV outlier variants

Finds variants in the ClusterBatch VCFs that are enriched for outlier
samples. The proportion of outlier samples in each variant is used to determine
outlier variants.
"""

import argparse
import sys
from pathlib import Path

import duckdb


# Flags to indicate the type of the outlier samples database being used.
# Outlier samples are from a custom list.
OUTLIER_SAMPLES_FROM_LIST = 0
# Outlier samples are from SV counts.
OUTLIER_SAMPLES_FROM_FILTERS = 1


def find_outliers_from_filter(
    con: duckdb.DuckDBPyConnection, filter_id: int, min_prop: float
):
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
        "SELECT DISTINCT l.vid"
        " FROM jrc_db.jrc_clusters l"
        f" JOIN var_db.outliers_{filter_id} r ON (l.member = r.vid);"
    )

    return [x[0] for x in con.sql(sql).fetchall()]


def find_outliers_from_list(
    con: duckdb.DuckDBPyConnection, i: int, svtype: str, min_prop: float
):
    sql = (
        f"CREATE OR REPLACE TABLE var_db.outliers_{i}"
        " AS SELECT vid FROM ("
        "SELECT l.vid AS vid,"
        " count(*) AS n_samples,"
        " count(*) FILTER (r.sample IS NOT NULL) AS n_outliers"
        " FROM ("
        "SELECT vid, sample"
        " FROM var_db.variants"
        " WHERE svtype = ?) l"
        " LEFT JOIN (SELECT sample FROM outlier_samples WHERE svtype = ?) r"
        " USING(sample)"
        " GROUP BY l.vid)"
        " WHERE n_outliers / n_samples >= ?;"
    )
    con.execute(sql, [svtype, svtype, min_prop])

    sql = (
        "SELECT DISTINCT l.vid"
        " FROM jrc_db.jrc_clusters l"
        f" JOIN var_db.outliers_{i} r ON (l.member = r.vid);"
    )

    return [x[0] for x in con.sql(sql).fetchall()]


def find_outlier_variants(
    con: duckdb.DuckDBPyConnection, min_prop: float, ols_dbtype: int
):
    if ols_dbtype == OUTLIER_SAMPLES_FROM_LIST:
        svtypes = con.sql("SELECT DISTINCT svtype FROM outlier_samples;").fetchall()
        outliers = [
            find_outliers_from_list(con, i, sv[0], min_prop)
            for i, sv in enumerate(svtypes)
        ]

        return [x for sublist in outliers for x in sublist]
    else:
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
        wgd_outliers = [x[0] for x in con.sql(sql).fetchall()]

        return [x for sublist in count_outliers for x in sublist] + wgd_outliers


def detect_outlier_samples_db_type(con: duckdb.DuckDBPyConnection):
    """Guess the type of outlier samples database that is connected."""
    tables = con.sql("SHOW TABLES;").fetchall()
    if len(tables) == 0:
        raise ValueError("Outlier samples database has no tables")
    if len(tables) == 1 and tables[0][0] == "outlier_samples":
        return OUTLIER_SAMPLES_FROM_LIST

    return OUTLIER_SAMPLES_FROM_FILTERS


def determine_outlier_variants(
    samples_db: Path, jrc_cluster_db: Path, variants_db: Path, min_prop: float
):
    with duckdb.connect(samples_db) as con:
        ols_dbtype = detect_outlier_samples_db_type(con)
        con.sql(f"ATTACH '{jrc_clusters_db}' AS jrc_db;")
        con.sql(f"ATTACH '{variants_db}' AS var_db;")
        outliers = find_outlier_variants(con, min_prop, ols_dbtype)

    for vid in outliers:
        print(vid)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Determine outlier variants in GATK-SV callset"
    )
    parser.add_argument(
        "outlier_samples_db",
        metavar="OUTLIER_SAMPLES_DB",
        help="Path to the outlier samples DuckDB database",
        type=Path,
    )
    parser.add_argument(
        "jrc_clusters_db",
        metavar="JOIN_RAW_CALLS_CLUSTERS_DB",
        help="Path to the JoinRawCalls clusters DuckDB database",
        type=Path,
    )
    parse.add_argument(
        "variants_db",
        metavar="VARIANTS_DB",
        help="Path to the ClusterBatch DuckDB database",
        type=Path,
    )
    parse.add_argument(
        "min_prop",
        metavar="MIN_PROP",
        help="Minimum proportion of outlier samples a variant needs to be considered an outlier",
        type=float,
    )

    retval = 0

    if not args.outlier_samples_db.is_file():
        raise FileNotFoundError("Outlier samples database not found")
    if not args.jrc_clusters_db.is_file():
        raise FileNotFoundError("JoinRawCalls clusters database not found")
    if not args.variants_db.is_file():
        raise FileNotFoundError("Variants database not found")
    if args.min_prop < 0 or args.min_prop > 1:
        raise ValueError("Min outlier sample proportion must be [0, 1]")

    determine_outlier_variants(
        args.outlier_samples_db, args.jrc_clusters_db, args.variants_db, args.min_prop
    )

    return retval


if __name__ == "__main__":
    raise SystemExit(main())
