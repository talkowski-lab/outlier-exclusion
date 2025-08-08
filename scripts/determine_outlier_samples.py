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

import pandas as pd


def find_sv_count_outliers(
    sv_counts: pd.DataFrame, filters: pd.DataFrame, iqr_mult: float
) -> pd.DataFrame:
    """Find outlier samples using SV counts."""
    filter_outliers = []
    for _, row in filters.iterrows():
        counts = sv_counts[
            (sv_counts["svtype"] == row["svtype"])
            & (sv_counts["svlen"] >= row["min_svlen"])
            & (sv_counts["svlen"] <= row["max_svlen"])
        ]
        if counts.empty:
            continue

        quantiles = counts["count"].quantile([0.25, 0.75])
        iqr = quantiles[0.75] - quantiles[0.25]
        median = counts["count"].median()
        lower_bound = median - iqr * iqr_mult
        upper_bound = median + iqr * iqr_mult
        
        outliers = counts[
            (counts["count"] < lower_bound) | (counts["count"] > upper_bound)
        ]
        outliers = outliers.assign(
            min_svlen=row["min_svlen"], max_svlen=row["max_svlen"]
        )
        filter_outliers.append(outliers)

    if not filter_outliers:
        return pd.DataFrame(
            columns=["sample", "count", "svtype", "min_svlen", "max_svlen"]
        )
    return pd.concat(filter_outliers)


def find_wgd_outliers(
    wgd_scores: pd.DataFrame, min_wgd: float, max_wgd: float
) -> pd.DataFrame:
    """Find outlier samples using WGD scores."""
    return wgd_scores[
        (wgd_scores["score"] < min_wgd) | (wgd_scores["score"] > max_wgd)
    ]


def determine_outliers(
    sv_counts_tsv: Path,
    sv_filters_tsv: Path,
    iqr_mult: float,
    sv_count_outliers_path: Path,
    wgd_outlier_samples_path: Path,
    wgd_scores: Path | None = None,
    min_wgd: float | None = None,
    max_wgd: float | None = None,
):
    sv_counts = pd.read_csv(sv_counts_tsv, sep="\t")
    sv_filters = pd.read_csv(sv_filters_tsv, sep="\t")
    sv_count_outliers = find_sv_count_outliers(sv_counts, sv_filters, iqr_mult)
    sv_count_outliers.to_csv(sv_count_outliers_path, sep="\t", index=False)
    
    if wgd_scores:
        if not wgd_scores.is_file():
            raise FileNotFoundError("WGD scores file must exist if given")
        if min_wgd is None or max_wgd is None:
            raise ValueError(
                "Min and max WGD scores must be given if WGD scores are given"
            )
        if min_wgd > max_wgd:
            raise ValueError("Min WGD score must be <= max WGD score")

        wgd_scores_df = pd.read_csv(wgd_scores, sep="\t", names=["sample", "score"])
        wgd_outliers = find_wgd_outliers(wgd_scores_df, min_wgd, max_wgd)
        wgd_outliers.to_csv(wgd_outlier_samples_path, sep="\t", index=False, header=False)
    else:
        # create empty file
        wgd_outlier_samples_path.touch()


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Determine outlier samples in GATK-SV callset"
    )
    parser.add_argument(
        "sv_counts_tsv",
        metavar="SV_COUNTS_TSV",
        help="Path to the input SV counts TSV",
        type=Path,
    )
    parser.add_argument(
        "sv_filters_tsv",
        metavar="SV_FILTERS_TSV",
        help="Path to the input SV filters TSV",
        type=Path,
    )
    parser.add_argument(
        "sv_count_outlier_samples_tsv",
        metavar="SV_COUNT_OUTLIER_SAMPLES_TSV",
        help="Output path for SV count outlier samples",
        type=Path,
    )
    parser.add_argument(
        "wgd_outlier_samples_tsv",
        metavar="WGD_OUTLIER_SAMPLES_TSV",
        help="Output path for WGD outlier samples",
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

    if not args.sv_counts_tsv.is_file():
        raise FileNotFoundError("Counts TSV must exist")
    if not args.sv_filters_tsv.is_file():
        raise FileNotFoundError("Filters TSV must exist")
    if args.iqr_mult < 0:
        raise ValueError("IQR multiplier must be greater than or equal to 0")

    determine_outliers(
        args.sv_counts_tsv,
        args.sv_filters_tsv,
        args.iqr_mult,
        args.sv_count_outlier_samples_tsv,
        args.wgd_outlier_samples_tsv,
        args.wgd_scores,
        args.min_wgd,
        args.max_wgd,
    )

    return retval


if __name__ == "__main__":
    raise SystemExit(main())
