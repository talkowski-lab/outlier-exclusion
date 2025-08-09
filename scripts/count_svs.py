"""Count SVs per sample

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

import argparse
from pathlib import Path
from collections.abc import Sequence

import pandas as pd


def validate_filters(filters: pd.DataFrame):
    if (filters["min_svlen"] < 0).any():
        raise ValueError("Min SV length must be >= 0")
    if (filters["min_svlen"] > filters["max_svlen"]).any():
        raise ValueError("Min SV length must be <= max SV length")


def count_svs(svs: pd.DataFrame, filters: pd.DataFrame) -> pd.DataFrame:
    """Count SVs per sample for each filter."""
    all_counts = []
    for _, row in filters.iterrows():
        filtered_svs = svs[
            (svs["svtype"] == row["svtype"])
            & (svs["svlen"] >= row["min_svlen"])
            & (svs["svlen"] <= row["max_svlen"])
        ]
        if filtered_svs.empty:
            continue
        counts = filtered_svs.groupby("sample").size().reset_index(name="count")
        counts = counts.assign(
            svtype=row["svtype"],
            min_svlen=row["min_svlen"],
            max_svlen=row["max_svlen"],
        )
        all_counts.append(counts)
    if not all_counts:
        return pd.DataFrame(
            columns=["sample", "count", "svtype", "min_svlen", "max_svlen"]
        )
    return pd.concat(all_counts, ignore_index=True)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Count SVs per sample")
    parser.add_argument(
        "svs_tsv",
        metavar="SVS_TSV",
        type=Path,
        help="Path to the SVs TSV file",
    )
    parser.add_argument(
        "filters_tsv",
        metavar="FILTERS_TSV",
        type=Path,
        help="Path to the filters TSV file",
    )
    parser.add_argument(
        "output_tsv",
        metavar="OUTPUT_TSV",
        type=Path,
        help="Path to the output TSV file",
    )
    args = parser.parse_args(argv)

    retval = 0

    if not args.svs_tsv.is_file():
        raise FileNotFoundError("SVs TSV file must exist")
    if not args.filters_tsv.is_file():
        raise FileNotFoundError("Filters TSV file must exist")

    svs = pd.read_csv(args.svs_tsv, sep="\t")
    filters = pd.read_csv(args.filters_tsv, sep="\t")

    validate_filters(filters)
    counts_df = count_svs(svs, filters)
    counts_df.to_csv(args.output_tsv, sep="\t", index=False)

    return retval


if __name__ == "__main__":
    raise SystemExit(main())
