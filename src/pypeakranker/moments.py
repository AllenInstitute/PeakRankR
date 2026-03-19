#!/usr/bin/env python3
"""
moments.py (pypeakranker)

Compute skewness, kurtosis, and bimodality coefficient per peak across multiple BigWig files.

Package API:
- add_moments(table_tsv, bigwig_files, out_tsv): read an existing feature table (chr/start/end),
  append per-bigwig metrics as new columns, and write out.

Standalone peaks -> moments table:
python -m pypeakranker.moments --peaks ... --bigwig-dir ... --output ...

Author: Saroja Somasundaram
Refactor: 2026-02 (table mode + peaks mode, CLI args, no hard-coded paths)
"""

from __future__ import annotations

import argparse
import os
from glob import glob
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
import pyBigWig
from scipy.stats import kurtosis, skew


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

    # add a stable peak_id for joining/pivoting in standalone mode
    df["peak_id"] = [f"peak_{i}" for i in range(len(df))]
    return df


def resolve_bigwigs(bigwig_dir: Optional[str], bigwig_files: Optional[List[str]]) -> List[str]:
    if (bigwig_dir is None) == (bigwig_files is None):
        raise SystemExit("Provide exactly one of --bigwig-dir OR --bigwig-files")

    if bigwig_dir is not None:
        bw_files = sorted(glob(os.path.join(bigwig_dir, "*.bw")))
    else:
        bw_files = sorted(list(bigwig_files or []))

    if not bw_files:
        raise FileNotFoundError("No BigWig files found.")
    return bw_files


def get_per_base_signal(
    bw: pyBigWig.pyBigWig,
    chrom: str,
    start: int,
    end: int,
    allow_missing_chroms: bool,
) -> np.ndarray:
    try:
        values = bw.values(chrom, int(start), int(end), numpy=True)
        return np.nan_to_num(values, nan=0.0)
    except RuntimeError:
        if allow_missing_chroms:
            return np.array([], dtype=float)
        raise
    except Exception:
        return np.array([], dtype=float)


def bimodality_coefficient(vals: np.ndarray) -> float:
    n = len(vals)
    if n < 4:
        return float("nan")
    s = skew(vals)
    k = kurtosis(vals, fisher=False)
    denom = k + (3 * (n - 1) ** 2) / ((n - 2) * (n - 3))
    if denom == 0:
        return float("nan")
    return float((s**2 + 1) / denom)


def _metrics_for_interval(vals: np.ndarray) -> Tuple[float, float, float]:
    if vals.size == 0:
        return (np.nan, np.nan, np.nan)
    if np.all(vals == 0) or np.all(vals == vals[0]):
        return (0.0, 0.0, np.nan)
    return (float(skew(vals)), float(kurtosis(vals, fisher=False)), bimodality_coefficient(vals))


# -------------------------
# Package API
# -------------------------

def add_moments(
    table_tsv: str,
    bigwig_files: List[str],
    out_tsv: str,
    allow_missing_chroms: bool = False,
    quiet: bool = False,
    prefix: str = "",
) -> None:
    """
    Read an existing feature table TSV (must contain chr/start/end),
    append per-bigwig skewness/kurtosis/bimodality columns, and write out.
    """
    base = pd.read_csv(table_tsv, sep="\t")
    needed = {"chr", "start", "end"}
    if not needed.issubset(base.columns):
        raise ValueError(f"Table must contain columns {needed}")

    if base.duplicated(["chr", "start", "end"]).any():
        raise ValueError("Table has duplicate (chr,start,end) rows; cannot merge safely.")

    # create columns in a deterministic order
    out = base.copy()

    for bw_path in bigwig_files:
        celltype = os.path.basename(bw_path).replace(".bw", "")
        col_skew = f"{prefix}{celltype}_skewness"
        col_kurt = f"{prefix}{celltype}_kurtosis"
        col_bim = f"{prefix}{celltype}_bimodality"

        log(f"Processing {celltype}", quiet)

        skew_vals: List[float] = []
        kurt_vals: List[float] = []
        bim_vals: List[float] = []

        with pyBigWig.open(bw_path) as bw:
            for _, row in out.iterrows():
                chrom = str(row["chr"])
                start = int(row["start"])
                end = int(row["end"])

                vals = get_per_base_signal(
                    bw=bw,
                    chrom=chrom,
                    start=start,
                    end=end,
                    allow_missing_chroms=allow_missing_chroms,
                )

                s, k, b = _metrics_for_interval(vals)
                skew_vals.append(s)
                kurt_vals.append(k)
                bim_vals.append(b)

        out[col_skew] = skew_vals
        out[col_kurt] = kurt_vals
        out[col_bim] = bim_vals

    ensure_parent_dir(out_tsv)
    out.to_csv(out_tsv, sep="\t", index=False)
    log(f"Wrote table with moments: {out_tsv}", quiet)


# -------------------------
# Standalone script CLI
# -------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Compute skewness, kurtosis, and bimodality per peak across BigWig files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Standalone peaks mode
    p.add_argument("--peaks", help="Peaks BED/TSV (headerless). Must have >=3 columns: chr, start, end.")

    # Table mode
    p.add_argument("--table", help="Existing feature table TSV (must have chr/start/end).")

    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--bigwig-dir", help="Directory containing *.bw BigWig files.")
    g.add_argument("--bigwig-files", nargs="+", help="Explicit list of BigWig files.")

    p.add_argument("--output", required=True, help="Output TSV path.")
    p.add_argument("--allow-missing-chroms", action="store_true")
    p.add_argument("--prefix", default="", help="Optional prefix for appended columns.")
    p.add_argument("--quiet", action="store_true")
    return p


def main() -> None:
    args = build_parser().parse_args()

    if (args.peaks is None) == (args.table is None):
        raise SystemExit("Provide exactly one of --peaks OR --table")

    bw_files = resolve_bigwigs(args.bigwig_dir, args.bigwig_files)
    log(f"Found {len(bw_files)} BigWig files.", args.quiet)

    if args.peaks:
        peaks_df = load_peaks(args.peaks, quiet=args.quiet)

        # Write peaks as a "table" then reuse add_moments logic for consistent output
        with tempfile_named_tsv(peaks_df) as tmp_table:
            add_moments(
                table_tsv=tmp_table,
                bigwig_files=bw_files,
                out_tsv=args.output,
                allow_missing_chroms=args.allow_missing_chroms,
                quiet=args.quiet,
                prefix=args.prefix,
            )
    else:
        add_moments(
            table_tsv=args.table,
            bigwig_files=bw_files,
            out_tsv=args.output,
            allow_missing_chroms=args.allow_missing_chroms,
            quiet=args.quiet,
            prefix=args.prefix,
        )


# helper context manager to avoid importing tempfile at top unless needed
class tempfile_named_tsv:
    def __init__(self, df: pd.DataFrame):
        self.df = df
        self.path = None

    def __enter__(self) -> str:
        import tempfile
        fd, path = tempfile.mkstemp(suffix=".tsv")
        os.close(fd)
        self.df.to_csv(path, sep="\t", index=False)
        self.path = path
        return path

    def __exit__(self, exc_type, exc, tb) -> None:
        try:
            if self.path and os.path.exists(self.path):
                os.remove(self.path)
        except Exception:
            pass


if __name__ == "__main__":
    main()