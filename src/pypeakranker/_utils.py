"""
_utils.py (pypeakranker)

Shared helpers used across submodules.
Not part of the public API.
"""

from __future__ import annotations

import os

import pandas as pd


def log(msg: str, quiet: bool) -> None:
    if not quiet:
        print(msg, flush=True)


def ensure_parent_dir(path: str) -> None:
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


def load_peaks(peaks_path: str, quiet: bool = False) -> pd.DataFrame:
    """Read a headerless peaks BED/TSV and name columns: chr, start, end, col4.."""
    log(f"Reading peaks from: {peaks_path}", quiet)
    df = pd.read_csv(peaks_path, sep="\t", header=None, comment="#")
    if df.shape[1] < 3:
        raise ValueError("Peaks file must have at least 3 columns: chr, start, end")

    n = df.shape[1]
    cols = ["chr", "start", "end"] + [f"col{i}" for i in range(4, n + 1)]
    df.columns = cols

    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    before = len(df)
    df = df.dropna(subset=["start", "end"]).copy()
    dropped = before - len(df)
    if dropped:
        log(f"Warning: Dropped {dropped} rows with non-numeric start/end.", quiet)

    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)

    bad = df["end"] <= df["start"]
    if bad.any():
        log(f"Warning: Dropping {bad.sum()} rows with end<=start.", quiet)
        df = df.loc[~bad].copy()

    dup = df.duplicated(subset=["chr", "start", "end"])
    if dup.any():
        log(f"Warning: Found {dup.sum()} duplicated peaks; removing duplicates.", quiet)
        df = df.drop_duplicates(subset=["chr", "start", "end"]).copy()

    return df
