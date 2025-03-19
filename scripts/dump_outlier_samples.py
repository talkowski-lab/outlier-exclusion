"""Dump outlier samples to text files

Dump the outlier sample tables in the DuckDB database to files. This is not
necessary for the pipeline, but can be helpful for users.
"""

import argparse
from pathlib import Path

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


def dump_sv_counts_outliers(db, outdir):
    with duckdb.connect(db) as con:
        filter_ids = con.sql("SELECT id FROM sv_filters;").fetchall()
        for i in filter_ids:
            outfile = Path(outdir, f"outliers_{i[0]}.tsv")
            dump_outliers(con, i[0], outfile)


def dump_wgd_outliers(db, path):
    # This seems like a SQL injection vunerability
    with duckdb.connect(db) as con:
        sql = f"COPY (SELECT sample, score FROM wgd_outliers LEFT JOIN wgd_scores USING(sample)) TO '{path}' (DELIMITER '\\t');"
        con.sql(sql)


def main(args):
    counts_db = Path(args.sv_counts_db)
    wgd_outliers = Path(args.wgd_outliers)
    outdir = Path(args.outdir)
    if not counts_db.is_file():
        raise FileNotFoundError("Counts database not found")
    if wgd_outliers.is_file():
        raise ValueError("WGD outliers file exists")
    if outdir.is_dir():
        raise ValueError("Output directory exists")
    outdir.mkdir()
    dump(counts_db, outdir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Dump outliers sample from DuckDB file"
    )
    parser.add_argument(
        "sv_counts_db",
        metavar="SV_COUNTS_DB",
        help="Path to the input SV counts DuckDB datatbase",
    )
    parse.add_argument(
        "wgd_outliers",
        metavar="WGD_OUTLIERS",
        help="Where to write the table of WGD outliers",
    )
    parser.add_argument("outdir", metavar="OUTDIR", help="Path to the output directory")

    main(parser.parse_args())
