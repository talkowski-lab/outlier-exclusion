"""Determine GATK-SV outlier samples

There are two types of outlier samples that can be defined: one based on SV
counts per sample and another based on WGD scores.  Outlier samples from SV
counts are determined per filter. For each SV type/size range, the SV counts
per sample are tabulated to compute an inter-quartile range which is then
multipled by a scaling factor to define a count limit. All samples with SV
counts exceeding the count limit are considered outlier samples. WGD outliers
are determined using fixed user-provided thresholds. All samples with scores
less than the given min or more than the given max are considered outliers.
"""

import argparse
from pathlib import Path
from collections.abc import Sequence

import duckdb


def find_filter_outliers(
    con: duckdb.DuckDBPyConnection, filter_id: int, iqr_mult: float
):
    sql = (
        "SELECT quant[2], quant[3] - quant[1]"
        " FROM"
        " (SELECT quantile_cont(count, [0.25, 0.5, 0.75]) AS quant"
        f" FROM sv_counts_{filter_id});"
    )
    results = con.sql(sql).fetchall()[0]
    median = results[0]
    iqr = results[1]
    sql = (
        f"CREATE OR REPLACE TABLE outliers_{filter_id}"
        " AS"
        " SELECT sample"
        f" FROM sv_counts_{filter_id}"
        " WHERE count < $1 - $2 * $3 OR count > $1 + $2 * $3"
    )
    con.execute(sql, [median, iqr, iqr_mult])


def find_sv_count_outliers(con: duckdb.DuckDBPyConnection, iqr_mult: float):
    filter_ids = con.sql("SELECT id FROM sv_filters;").fetchall()
    for i in filter_ids:
        find_filter_outliers(con, i[0], iqr_mult)


def find_wgd_outliers(con: duckdb.DuckDBPyConnection, min_wgd: float, max_wgd: float):
    sql = (
        "CREATE OR REPLACE TABLE wgd_outliers"
        " AS SELECT sample"
        " FROM wgd_scores"
        " WHERE score < ? OR score > ?"
    )
    con.execute(sql, [min_wgd, max_wgd])


def make_wgd_table(con: duckdb.DuckDBPyConnection):
    con.sql("CREATE OR REPLACE TABLE wgd_scores (sample VARCHAR, score FLOAT);")


def load_wgd_table(con: duckdb.DuckDBPyConnection, scores_path: Path):
    # This seems like a SQL injection vunerability
    con.sql(f"COPY wgd_scores FROM '{scores_path}' (DELIMITER '\\t', HEADER false);")


def determine_outliers(
    db: Path,
    iqr_mult: float,
    wgd_scores: Path | None = None,
    min_wgd: float | None = None,
    max_wgd: float | None = None,
):
    # TODO Might be sensible to do this in a transaction and roll back on error
    with duckdb.connect(db) as con:
        find_sv_count_outliers(con, iqr_mult)
        make_wgd_table(con)
        if wgd_scores is not None:
            if not wgd_scores.is_file():
                raise FileNotFoundError("WGD scores file must exist if given")
            else:
                load_wgd_table(con, wgd_scores)
            if min_wgd is None or max_wgd is None:
                raise ValueError(
                    "Min and max WGD scores must be given if WGD scores are given"
                )
            find_wgd_outliers(con, min_wgd, max_wgd)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Determine outlier samples in GATK-SV callset"
    )
    parser.add_argument(
        "sv_counts_db",
        metavar="SV_COUNTS_DB",
        help="Path to the input SV counts DuckDB database",
        type=Path,
    )
    parser.add_argument(
        "iqr_mult",
        metavar="IQR_MULTIPLIER",
        help="SVs per genome IQR multiplier",
        type=float,
    )
    parser.add_argument(
        "--wgd-scores",
        dest="wgd_scores",
        metavar="WGD_SCORES",
        help="Path to the sample WGD scores",
        type=Path,
    )
    parser.add_argument(
        "--min-wgd",
        dest="min_wgd",
        metavar="MIN_WGD",
        help="Minimum WGD score. Must be given if --wgd-scores is given",
        type=float,
    )
    parser.add_argument(
        "--max-wgd",
        dest="max_wgd",
        metavar="MAX_WGD",
        help="Maximum WGD score. Must be given if --wgd-scores is given",
        type=float,
    )
    args = parser.parse_args(argv)

    retval = 0

    if not args.sv_counts_db.is_file():
        raise FileNotFoundError("Counts database must exist")
    if args.iqr_mult < 0:
        raise ValueError("IQR multiplier must be greater than or equal to 0")
    if args.min_wgd > args.max_wgd:
        raise ValueError("Min WGD score must be >= max wgd score")

    determine_outliers(
        args.sv_counts_db, args.iqr_mult, args.wgd_scores, args.min_wgd, args.max_wgd
    )

    return retval


if __name__ == "__main__":
    raise SystemExit(main())
