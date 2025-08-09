"""Determine GATK-SV outlier variants

Finds variants in the ClusterBatch VCFs that are enriched for outlier
samples. The proportion of outlier samples in each variant is used to determine
outlier variants.
"""

import argparse
from pathlib import Path
from collections.abc import Sequence

import pandas as pd


def find_outlier_variants(
    variants: pd.DataFrame,
    clusters: pd.DataFrame,
    outlier_samples: pd.DataFrame,
    min_prop: float,
) -> pd.DataFrame:
    """Find outlier variants."""
    # Merge with clusters to get the cluster ID for each variant
    variants_with_clusters = pd.merge(variants, clusters, left_on="vid", right_on="member")

    # Merge with outlier samples to identify which carriers are outliers
    merged = pd.merge(
        variants_with_clusters,
        outlier_samples,
        on=["sample", "svtype"],
        how="left",
        indicator=True,
    )
    merged["is_outlier"] = merged["_merge"] == "both"

    # Count total samples and outlier samples per variant
    summary = (
        merged.groupby("vid_y")
        .agg(n_samples=("sample", "size"), n_outliers=("is_outlier", "sum"))
        .reset_index()
    )
    summary["outlier_prop"] = summary["n_outliers"] / summary["n_samples"]

    # Filter for outlier variants
    outlier_variants = summary[summary["outlier_prop"] >= min_prop]
    return outlier_variants


def determine_outlier_variants(
    variants_tsv: Path,
    clusters_tsv: Path,
    outlier_samples_tsv: Path,
    output: Path,
    min_prop: float,
):
    variants = pd.read_csv(variants_tsv, sep="\t",
                             names=["vid", "svtype", "svlen", "sample"])
    clusters = pd.read_csv(clusters_tsv, sep="\t", names=["vid", "member"])
    outlier_samples = pd.read_csv(outlier_samples_tsv, sep="\t")

    # The logic is simplified to handle one type of outlier source.
    # The original script had flags for different outlier sample sources.
    # This implementation assumes a single, consolidated list of outlier samples.
    outlier_variants_df = find_outlier_variants(
        variants, clusters, outlier_samples, min_prop
    )

    with open(output, "w") as f:
        for vid in set(outlier_variants_df["vid_y"]):
            f.write(vid + "\n")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Determine outlier variants in GATK-SV callset"
    )
    parser.add_argument(
        "variants_tsv",
        metavar="VARIANTS_TSV",
        help="Path to the variants TSV file",
        type=Path,
    )
    parser.add_argument(
        "jrc_clusters_tsv",
        metavar="JOIN_RAW_CALLS_CLUSTERS_TSV",
        help="Path to the JoinRawCalls clusters TSV file",
        type=Path,
    )
    parser.add_argument(
        "outlier_samples_tsv",
        metavar="OUTLIER_SAMPLES_TSV",
        help="Path to the outlier samples TSV file",
        type=Path,
    )
    parser.add_argument(
        "output",
        metavar="OUTPUT",
        help="Where to write the outlier variants",
        type=Path,
    )
    parser.add_argument(
        "min_prop",
        metavar="MIN_PROP",
        help="Minimum proportion of outlier samples a variant needs to be considered an outlier",
        type=float,
    )
    args = parser.parse_args(argv)

    retval = 0

    if not args.variants_tsv.is_file():
        raise FileNotFoundError("Variants TSV not found")
    if not args.jrc_clusters_tsv.is_file():
        raise FileNotFoundError("JoinRawCalls clusters TSV not found")
    if not args.outlier_samples_tsv.is_file():
        raise FileNotFoundError("Outlier samples TSV not found")
    if args.min_prop < 0 or args.min_prop > 1:
        raise ValueError("Min outlier sample proportion must be [0, 1]")
    if args.output.is_file():
        raise ValueError("Output file must not exist")

    determine_outlier_variants(
        args.variants_tsv,
        args.jrc_clusters_tsv,
        args.outlier_samples_tsv,
        args.output,
        args.min_prop,
    )

    return retval


if __name__ == "__main__":
    raise SystemExit(main())
