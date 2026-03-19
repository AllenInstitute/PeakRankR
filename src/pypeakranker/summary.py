#!/usr/bin/env python3
"""
summary.py (pypeakranker)

Library + CLI-compatible BigWig summarization.

- init_table(peaks_path, out_tsv): create a feature table from peaks (chr/start/end + col4..)
- add_signal(table_tsv, bigwig_files, out_tsv, ...): append bigWig summary columns to an existing table

As a standalone script:
python -m pypeakranker.summary --peaks ... --bigwig-files ... --output ...
"""

from __future__ import annotations

import argparse
import os
from glob import glob
from typing import List, Optional

import numpy as np
import pandas as pd
import pyBigWig


# -------------------------
# Shared helpers
# -------------------------

def log(msg: str, quiet: bool) -> None:
    if not quiet:
        print(msg, flush=True)


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


def ensure_parent_dir(path: str) -> None:
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


def sample_name(path: str, mode: str = "stem") -> str:
    base = os.path.basename(path)
    if mode == "filename":
        return base
    for ext in (".bw", ".bigWig", ".bigwig"):
        if base.endswith(ext):
            return base[: -len(ext)]
    return os.path.splitext(base)[0]


def summarize_values(values: np.ndarray, stat: str, keep_nans: bool) -> float:
    if values is None:
        return 0.0
    arr = np.asarray(values, dtype=float)
    if not keep_nans:
        arr = arr[~np.isnan(arr)]
    if arr.size == 0:
        return 0.0
    if stat == "sum":
        return float(arr.sum())
    if stat == "mean":
        return float(arr.mean())
    if stat == "max":
        return float(arr.max())
    raise ValueError(f"Unsupported stat: {stat}")


def summarize_peak(
    bw: pyBigWig.pyBigWig,
    chrom: str,
    start: int,
    end: int,
    stat: str,
    keep_nans: bool,
    allow_missing_chroms: bool,
) -> float:
    try:
        if allow_missing_chroms and chrom not in bw.chroms():
            return 0.0
        vals = bw.values(chrom, start, end, numpy=True)
        return summarize_values(vals, stat=stat, keep_nans=keep_nans)
    except RuntimeError:
        if allow_missing_chroms:
            return 0.0
        raise
    except Exception:
        # permissive: keep pipeline moving
        return 0.0


# -------------------------
# Package API
# -------------------------

def init_table(peaks_path: str, out_tsv: str, quiet: bool = False) -> None:
    """Initialize a feature table TSV from peaks."""
    df = load_peaks(peaks_path, quiet=quiet)
    ensure_parent_dir(out_tsv)
    df.to_csv(out_tsv, sep="\t", index=False)
    log(f"Wrote init table: {out_tsv}", quiet)


def add_signal(
    table_tsv: str,
    bigwig_files: List[str],
    out_tsv: str,
    stat: str = "sum",
    suffix: str = "summary",
    sample_name_mode: str = "stem",
    keep_nans: bool = False,
    allow_missing_chroms: bool = False,
    quiet: bool = False,
) -> None:
    """
    Read an existing feature table TSV (must contain chr/start/end),
    append one column per bigWig, and write out.
    """
    base = pd.read_csv(table_tsv, sep="\t")
    needed = {"chr", "start", "end"}
    if not needed.issubset(base.columns):
        raise ValueError(f"Table must contain columns {needed}")

    if base.duplicated(["chr", "start", "end"]).any():
        raise ValueError("Table has duplicate (chr,start,end) rows; cannot merge safely.")

    out = base.copy()

    for bw_path in bigwig_files:
        samp = sample_name(bw_path, mode=sample_name_mode)
        col = f"{samp}_{suffix}"
        log(f"Processing {samp} -> {col}", quiet)

        with pyBigWig.open(bw_path) as bw:
            out[col] = base.apply(
                lambda r: summarize_peak(
                    bw=bw,
                    chrom=str(r["chr"]),
                    start=int(r["start"]),
                    end=int(r["end"]),
                    stat=stat,
                    keep_nans=keep_nans,
                    allow_missing_chroms=allow_missing_chroms,
                ),
                axis=1,
            )

    ensure_parent_dir(out_tsv)
    out.to_csv(out_tsv, sep="\t", index=False)
    log(f"Wrote table with signal columns: {out_tsv}", quiet)


# -------------------------
# Standalone script CLI
# -------------------------

def _list_bigwigs(bigwig_dir: Optional[str], bigwig_files: Optional[List[str]], pattern: str) -> List[str]:
    if bigwig_files:
        files = list(bigwig_files)
    else:
        assert bigwig_dir is not None
        patterns = [p.strip() for p in pattern.split(",") if p.strip()]
        files = []
        for pat in patterns:
            files.extend(glob(os.path.join(bigwig_dir, pat)))
    files = sorted(set(files))
    if not files:
        raise FileNotFoundError("No BigWig files found. Check --bigwig-dir/--bigwig-files and --pattern.")
    return files

# -------------------------
# Standalone script CLI
# -------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Summarize BigWig signal over peaks (sum/mean/max).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    p.add_argument("--peaks", help="Peaks BED/TSV (headerless).")
    p.add_argument("--table", help="Existing feature table TSV (must have chr/start/end).")

    bw_src = p.add_mutually_exclusive_group(required=True)
    bw_src.add_argument("--bigwig-dir", help="Directory containing BigWig files.")
    bw_src.add_argument("--bigwig-files", nargs="+", help="One or more BigWig file paths.")

    p.add_argument("--pattern", default="*.bw,*.bigWig,*.bigwig", help="Glob patterns for --bigwig-dir.")
    p.add_argument("--output", required=True, help="Output TSV path.")
    p.add_argument("--stat", choices=["sum", "mean", "max"], default="sum")
    p.add_argument("--sample-name-mode", choices=["stem", "filename"], default="stem")
    p.add_argument("--suffix", default="summary", help="Name columns <sample>_<suffix>.")
    p.add_argument("--keep-nans", action="store_true")
    p.add_argument("--allow-missing-chroms", action="store_true")
    p.add_argument("--quiet", action="store_true")
    return p


def main() -> None:
    args = build_parser().parse_args()

    if (args.peaks is None) == (args.table is None):
        raise SystemExit("Provide exactly one of --peaks OR --table")

    bw_files = _list_bigwigs(args.bigwig_dir, args.bigwig_files, args.pattern)

    if args.peaks:
        # peaks -> init table -> add_signal
        import tempfile
        fd, tmp = tempfile.mkstemp(suffix=".tsv")
        os.close(fd)
        try:
            init_table(args.peaks, tmp, quiet=args.quiet)
            add_signal(
                table_tsv=tmp,
                bigwig_files=bw_files,
                out_tsv=args.output,
                stat=args.stat,
                suffix=args.suffix,
                sample_name_mode=args.sample_name_mode,
                keep_nans=args.keep_nans,
                allow_missing_chroms=args.allow_missing_chroms,
                quiet=args.quiet,
            )
        finally:
            try:
                os.remove(tmp)
            except Exception:
                pass
    else:
        # table -> add_signal
        add_signal(
            table_tsv=args.table,
            bigwig_files=bw_files,
            out_tsv=args.output,
            stat=args.stat,
            suffix=args.suffix,
            sample_name_mode=args.sample_name_mode,
            keep_nans=args.keep_nans,
            allow_missing_chroms=args.allow_missing_chroms,
            quiet=args.quiet,
        )


if __name__ == "__main__":
    main()