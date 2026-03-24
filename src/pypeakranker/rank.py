#!/usr/bin/env python3
"""
rank.py (pypeakranker)

Rank peaks by ATAC-seq signal specificity.

For each peak, specificity is defined as:

    ratio = target_signal / mean(all_background_signals)

where target_signal and background_signal values are taken from named
signal columns in the feature table (produced by add-signal).
Ratios are min-max normalised to [0, 1] across all peaks.

A higher specificity score means the peak has stronger signal in the
target group relative to all other groups.

Package API:
    rank_by_specificity(table_tsv, target_cols, background_cols, out_tsv, ...)

CLI:
    pypeakranker rank-specificity \
        --table features.tsv \
        --target-cols Astrocytes-1_mean \
        --background-cols Astrocytes-1_mean Astrocytes-2_mean Neurons_mean \
        --out ranked.tsv
"""

from __future__ import annotations

import argparse
import os
from typing import List, Optional

import numpy as np
import pandas as pd
from pypeakranker._utils import log, ensure_parent_dir


def _minmax_normalise(series: pd.Series) -> pd.Series:
    """Min-max normalise a Series to [0, 1]. Returns 0 everywhere if constant."""
    lo = series.min()
    hi = series.max()
    if hi == lo:
        return pd.Series(np.zeros(len(series)), index=series.index)
    return (series - lo) / (hi - lo)


def rank_by_specificity(
    table_tsv: str,
    target_cols: List[str],
    background_cols: List[str],
    out_tsv: str,
    target_agg: str = "mean",
    specificity_col: str = "specificity_score",
    rank_col: str = "specificity_rank",
    ascending: bool = False,
    quiet: bool = False,
) -> pd.DataFrame:
    """
    Read an existing feature table, compute per-peak ATAC specificity,
    and write a ranked output table.

    Parameters
    ----------
    table_tsv : str
        Path to feature table TSV (must contain chr/start/end and signal columns).
    target_cols : list of str
        Column name(s) representing the target group's ATAC signal.
        If multiple columns are given, they are aggregated (mean by default).
    background_cols : list of str
        Column names representing ALL groups' ATAC signal (including target).
        Used to compute average background signal per peak.
    out_tsv : str
        Output TSV path.
    target_agg : str
        How to aggregate multiple target columns: "mean" (default) or "sum".
    specificity_col : str
        Name of the appended specificity score column.
    rank_col : str
        Name of the appended rank column.
    ascending : bool
        If True, rank 1 = lowest score. Default False (rank 1 = highest score).
    quiet : bool
        Suppress progress messages.

    Returns
    -------
    pd.DataFrame
        The output table with specificity_score and specificity_rank appended.
    """
    df = pd.read_csv(table_tsv, sep="\t")
    needed = {"chr", "start", "end"}
    if not needed.issubset(df.columns):
        raise ValueError(f"Table must contain columns {needed!r}")

    # Validate columns exist
    missing_target = [c for c in target_cols if c not in df.columns]
    if missing_target:
        raise ValueError(
            f"Target column(s) not found in table: {missing_target}\n"
            f"Available columns: {list(df.columns)}"
        )

    missing_bg = [c for c in background_cols if c not in df.columns]
    if missing_bg:
        raise ValueError(
            f"Background column(s) not found in table: {missing_bg}\n"
            f"Available columns: {list(df.columns)}"
        )

    if not target_cols:
        raise ValueError("target_cols must not be empty.")
    if not background_cols:
        raise ValueError("background_cols must not be empty.")

    log(f"Target columns    : {target_cols}", quiet)
    log(f"Background columns: {background_cols}", quiet)
    log(f"Peaks             : {len(df)}", quiet)

    # ── Compute target signal ──────────────────────────────────────────────
    if len(target_cols) == 1:
        target_signal = df[target_cols[0]].fillna(0.0)
    elif target_agg == "mean":
        target_signal = df[target_cols].fillna(0.0).mean(axis=1)
    elif target_agg == "sum":
        target_signal = df[target_cols].fillna(0.0).sum(axis=1)
    else:
        raise ValueError(f"target_agg must be 'mean' or 'sum', got {target_agg!r}")

    # ── Compute per-peak average background signal ─────────────────────────
    # background_cols includes target group(s) — average across all groups
    avg_bg = df[background_cols].fillna(0.0).mean(axis=1)

    # ── Specificity ratio ──────────────────────────────────────────────────
    # ratio = target / avg_bg; peaks with zero background get ratio 0
    ratio = target_signal.copy()
    nonzero = avg_bg > 0
    ratio[nonzero] = target_signal[nonzero] / avg_bg[nonzero]
    ratio[~nonzero] = 0.0

    # ── Min-max normalise ──────────────────────────────────────────────────
    spec_score = _minmax_normalise(ratio)

    # ── Rank (1 = best by default = highest specificity) ───────────────────
    spec_rank = spec_score.rank(method="min", ascending=ascending).astype(int)

    out = df.copy()
    out[specificity_col] = spec_score.values
    out[rank_col]        = spec_rank.values

    # Sort by rank
    out = out.sort_values(rank_col).reset_index(drop=True)

    ensure_parent_dir(out_tsv)
    out.to_csv(out_tsv, sep="\t", index=False)
    log(f"Wrote ranked table: {out_tsv}", quiet)
    log(f"Top peak: {out.iloc[0]['chr']}:{out.iloc[0]['start']}-{out.iloc[0]['end']}  "
        f"specificity={out.iloc[0][specificity_col]:.4f}", quiet)

    return out


# ─────────────────────────────────────────────────────────────────────────────
# Standalone CLI
# ─────────────────────────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Rank peaks by ATAC-seq specificity.\n\n"
            "Specificity is defined per peak as:\n"
            "  ratio = target_signal / mean(all_background_signals)\n"
            "Ratios are min-max normalised to [0, 1] across all peaks.\n"
            "Rank 1 = most specific to the target group."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--table", required=True,
                   help="Feature table TSV (must contain chr/start/end and signal columns).")
    p.add_argument("--target-cols", nargs="+", required=True,
                   help="Signal column(s) for the TARGET group (e.g. Astrocytes-1_mean).")
    p.add_argument("--background-cols", nargs="+", required=True,
                   help="Signal columns for ALL groups including target "
                        "(e.g. Astrocytes-1_mean Neurons_mean Oligo_mean).")
    p.add_argument("--out", required=True,
                   help="Output TSV path.")
    p.add_argument("--target-agg", default="mean", choices=["mean", "sum"],
                   help="How to aggregate multiple --target-cols. Default: mean.")
    p.add_argument("--specificity-col", default="specificity_score",
                   help="Name of the appended specificity score column.")
    p.add_argument("--rank-col", default="specificity_rank",
                   help="Name of the appended rank column.")
    p.add_argument("--quiet", action="store_true")
    return p


def main() -> None:
    args = build_parser().parse_args()
    rank_by_specificity(
        table_tsv       = args.table,
        target_cols     = args.target_cols,
        background_cols = args.background_cols,
        out_tsv         = args.out,
        target_agg      = args.target_agg,
        specificity_col = args.specificity_col,
        rank_col        = args.rank_col,
        quiet           = args.quiet,
    )


if __name__ == "__main__":
    main()