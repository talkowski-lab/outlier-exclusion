"""Dump outlier samples to text files

Dump the outlier sample tables in the DuckDB database to files. This is not
necessary for the pipeline, but can be helpful for users.
"""

import argparse
from pathlib import Path
from collections.abc import Sequence

import duckdb


def dump_filter_outliers(con, filter_id, outfile):
    sql = (
        "COPY (SELECT sample, count, svtype, min_svlen, max_svlen"
        " FROM ("
        f"SELECT sample, count, {filter_id} AS id"
        f" FROM outliers_{filter_id}"
        f" LEFT JOIN sv_counts_{filter_id}"
        " USING (sample))"
        " JOIN sv_filters USING (id))"
        f" TO '{outfile}'"
        " (DELIMITER '\t', HEADER false);"
    )
    con.sql(sql)


def dump_sv_counts_outliers(db: Path, outdir: Path):
    with duckdb.connect(db) as con:
        filter_ids = con.sql("SELECT id FROM sv_filters;").fetchall()
        for i in filter_ids:
            outfile = Path(outdir, f"outliers_{i[0]}.tsv")
            dump_filter_outliers(con, i[0], outfile)


def dump_wgd_outliers(db: Path, outfile: Path):
    # This seems like a SQL injection vunerability
    with duckdb.connect(db) as con:
        sql = f"COPY (SELECT sample, score FROM wgd_outliers LEFT JOIN wgd_scores USING(sample)) TO '{outfile}' (DELIMITER '\\t');"
        con.sql(sql)


def dump(counts_db: Path, wgd_outliers: Path, outdir: Path):
    outdir.mkdir()
    dump_sv_counts_outliers(counts_db, outdir)
    dump_wgd_outliers(counts_db, wgd_outliers)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Dump outliers sample from DuckDB database"
    )
    parser.add_argument(
        "counts_db",
        metavar="SV_COUNTS_DB",
        type=Path,
        help="Path to the input SV counts DuckDB datatbase",
    )
    parser.add_argument(
        "wgd_outliers",
        metavar="WGD_OUTLIERS",
        type=Path,
        help="Where to write the table of WGD outliers",
    )
    parser.add_argument(
        "outdir", metavar="OUTDIR", type=Path, help="Path to the output directory"
    )
    args = parser.parse_args(argv)

    retval = 0

    if not args.counts_db.is_file():
        raise FileNotFoundError("Counts database not found")
    if args.wgd_outliers.is_file():
        raise ValueError("WGD outliers file exists")
    if args.outdir.is_dir():
        raise ValueError("Output directory exists")

    dump(args.counts_db, args.wgd_outliers, args.outdir)

    return retval


if __name__ == "__main__":
    raise SystemExit(main())
