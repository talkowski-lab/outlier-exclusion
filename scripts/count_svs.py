"""Count SVs per sample in a DuckDB database

usage: python count_svs.py <counts_db> <sv_db>

Given a database of filters in <counts_db> and a database of SVs in <sv_db>,
count the number of SVs per sample in <sv_db> corresponding to the filters in
<counts_db>. Each filter defines a SV type and SV size range combination to
count.

==================
| schema <sv_db> |
==================
CREATE TABLE svs (
    vid VARCHAR,
    svtype VARCHAR,
    svlen INTEGER,
    sample VARCHAR
);

======================
| schema <counts_db> |
======================
CREATE TABLE sv_filters (
    svtype VARCHAR,
    min_svlen DOUBLE,
    max_svlen DOUBLE
);
CREATE SEQUENCE id_sequence START 1;
ALTER TABLE sv_filters ADD COLUMN id INTEGER DEFAULT nextval('id_sequence');


For each 'id' in the 'sv_filters' table, a table of the form
'sv_counts_{id}' will be written to <counts_db> containing the SV counts
per sample for the corresponding filter.

EXISTING TABLES WILL BE OVERWRITTEN!
"""

import sys
import argparse
from pathlib import Path
from collections.abc import Sequence

import duckdb


def validate_filters(con: duckdb.DuckDBPyConnection):
    sql = "SELECT svtype, min_svlen, max_svlen FROM sv_filters;"
    filters = con.sql(sql).fetchall()
    for f in filters:
        if f[1] < 0:
            raise ValueError("Min SV length must be >= 0")
        if f[1] > f[2]:
            raise ValueError("Min SV length must be <= max SV length")


def count_svs(con: duckdb.DuckDBPyConnection, filter_id: int):
    sql = """SELECT svtype, min_svlen, max_svlen
    FROM sv_filters
    WHERE id = ?;
    """
    filters = con.execute(sql, [filter_id]).fetchall()
    sql = (
        f"CREATE OR REPLACE TABLE sv_counts_{filter_id}"
        " AS SELECT sample, COUNT(*) AS count"
        " FROM sv_db.svs"
        " WHERE svtype = ? AND svlen >= ? AND svlen <= ?"
        " GROUP BY sample"
    )

    con.execute(sql, list(filters[0]))


def make_tables(counts_db: Path, sv_db: Path):
    with duckdb.connect(counts_db) as con:
        validate_filters(con)
        con.sql(f"ATTACH '{sv_db}' AS sv_db;")
        filter_ids = con.sql("SELECT id FROM sv_filters;").fetchall()
        for i in filter_ids:
            count_svs(con, i[0])


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Count SVs per sample")
    parser.add_argument(
        "counts_db",
        meta="COUNTS_DB",
        type=Path,
        help="Path to the SV counts DuckDB database",
    )
    parser.add_argument(
        "sv_db", meta="SV_DB", type=Path, help="Path to the SV DuckDB database"
    )
    args = parser.parse_args(argv)

    retval = 0

    if not args.counts_db.is_file():
        raise FileNotFoundError("Counts database must exist")
    if not args.sv_db.is_file():
        raise FileNotFoundError("SV database must exist")

    make_tables(args.counts_db, args.sv_db)

    return retval


if __name__ == "__main__":
    raise SystemExit(main())
