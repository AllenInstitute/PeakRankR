#!/usr/bin/env python3
"""
phylop.py (pypeakranker)

Compute mean PhyloP score per peak from a PhyloP bigWig track.

Modes:
- add_phylop(table_tsv, phylop_bw, out_tsv): read an existing feature table (chr/start/end),
  append phyloP_mean, write out (can overwrite same file).

Standalone peaks -> phyloP table:
python -m pypeakranker.phylop --peaks ... --phylop-bw ... --output ...

Optional liftOver:
If peaks/table are in a SOURCE assembly and the PhyloP bigWig is in a TARGET assembly,
provide --chain (sourceToTarget.over.chain[.gz]) and (optionally) --liftover-exe.

Example:
python -m pypeakranker.phylop \
  --table features.tsv \
  --phylop-bw hg38.phyloP100way.bw \
  --chain rheMac10ToHg38.over.chain.gz \
  --output features_with_phylop.tsv
"""

from __future__ import annotations

import argparse
import os
import subprocess
import tempfile
from typing import Optional, Tuple

import pandas as pd
import pyBigWig


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


def _run_liftover(
    peaks_df: pd.DataFrame,
    chain_file: str,
    liftover_exe: str,
    quiet: bool,
) -> pd.DataFrame:
    """
    Run UCSC liftOver on peaks_df and return a DataFrame with lifted coords.

    Input columns required: chr, start, end
    Output columns: chr_target, start_target, end_target
    Unmapped peaks will have NA in target columns.
    """
    for c in ["chr", "start", "end"]:
        if c not in peaks_df.columns:
            raise ValueError(f"liftOver requires column '{c}'")

    with tempfile.TemporaryDirectory() as tmp:
        in_bed = os.path.join(tmp, "in.bed")
        out_bed = os.path.join(tmp, "out.bed")
        un_bed = os.path.join(tmp, "unmapped.bed")

        # Use row index as name so we can merge back 1:1
        bed = peaks_df[["chr", "start", "end"]].copy()
        bed["name"] = peaks_df.index.astype(str)
        bed.to_csv(in_bed, sep="\t", header=False, index=False)

        cmd = [liftover_exe, in_bed, chain_file, out_bed, un_bed]
        log(f"Running liftOver: {' '.join(cmd)}", quiet)

        try:
            subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"liftOver failed.\nCMD: {' '.join(cmd)}\nSTDOUT:\n{e.stdout}\nSTDERR:\n{e.stderr}"
            )

        if not os.path.exists(out_bed) or os.path.getsize(out_bed) == 0:
            # nothing mapped
            out = peaks_df.copy()
            out["chr_target"] = pd.NA
            out["start_target"] = pd.NA
            out["end_target"] = pd.NA
            return out

        mapped = pd.read_csv(out_bed, sep="\t", header=None)
        mapped = mapped.iloc[:, :4]
        mapped.columns = ["chr_target", "start_target", "end_target", "name"]
        mapped["name"] = mapped["name"].astype(str)

        out = peaks_df.copy()
        out["__rowname__"] = out.index.astype(str)

        out = out.merge(
            mapped,
            left_on="__rowname__",
            right_on="name",
            how="left",
        ).drop(columns=["__rowname__", "name"])

        # ensure numeric start/end if present
        for c in ["start_target", "end_target"]:
            if c in out.columns:
                out[c] = pd.to_numeric(out[c], errors="coerce")

        return out


def _mean_bigwig(bw: pyBigWig.pyBigWig, chrom: str, start: int, end: int) -> float:
    try:
        vals = bw.values(chrom, start, end)
        vals = [v for v in vals if v is not None]
        return float(sum(vals) / len(vals)) if vals else 0.0
    except Exception:
        return 0.0


# -------------------------
# Package API
# -------------------------

def add_phylop(
    table_tsv: str,
    phylop_bw: str,
    out_tsv: str,
    chain_file: Optional[str] = None,
    liftover_exe: str = "liftOver",
    allow_missing_chroms: bool = False,
    max_len: int = 5000,
    quiet: bool = False,
    out_col: str = "phyloP_mean",
    keep_lifted_coords: bool = True,
) -> None:
    """
    Read an existing feature table TSV (must contain chr/start/end),
    append phyloP_mean, and write out.

    If chain_file is provided, liftOver peaks to the target assembly first,
    then query PhyloP on (chr_target, start_target, end_target).
    """
    base = pd.read_csv(table_tsv, sep="\t")
    needed = {"chr", "start", "end"}
    if not needed.issubset(base.columns):
        raise ValueError(f"Table must contain columns {needed}")

    if base.duplicated(["chr", "start", "end"]).any():
        raise ValueError("Table has duplicate (chr,start,end) rows; cannot merge safely.")

    work = base.copy()

    if chain_file:
        work = _run_liftover(work, chain_file=chain_file, liftover_exe=liftover_exe, quiet=quiet)
        chr_c, start_c, end_c = "chr_target", "start_target", "end_target"
    else:
        chr_c, start_c, end_c = "chr", "start", "end"

    log(f"Opening PhyloP bigWig: {phylop_bw}", quiet)
    bw = pyBigWig.open(phylop_bw)

    scores = []
    missing = 0
    too_long = 0

    for _, row in work.iterrows():
        chrom = row.get(chr_c)
        if pd.isna(chrom):
            missing += 1
            scores.append(0.0)
            continue

        start = int(row[start_c])
        end = int(row[end_c])

        if end <= start:
            scores.append(0.0)
            continue

        if (end - start) > max_len:
            too_long += 1
            scores.append(0.0)
            continue

        # allow_missing_chroms: if chrom not present in BW, return 0 instead of raising
        try:
            scores.append(_mean_bigwig(bw, str(chrom), start, end))
        except RuntimeError:
            missing += 1
            if allow_missing_chroms:
                scores.append(0.0)
            else:
                raise

    bw.close()

    if missing and not quiet:
        log(f"Warning: {missing} peaks missing target coordinates or chrom not in bigWig; filled as 0.", quiet)
    if too_long and not quiet:
        log(f"Warning: {too_long} peaks exceeded max_len={max_len}; filled as 0.", quiet)

    out = work.copy()
    out[out_col] = scores

    if chain_file and not keep_lifted_coords:
        out = out.drop(columns=["chr_target", "start_target", "end_target"], errors="ignore")

    ensure_parent_dir(out_tsv)
    out.to_csv(out_tsv, sep="\t", index=False)
    log(f"Wrote table with {out_col}: {out_tsv}", quiet)


# -------------------------
# Standalone script CLI
# -------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Compute mean PhyloP per peak from a PhyloP bigWig (optionally using liftOver).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Backward-compatible standalone mode (peaks -> phyloP table)
    p.add_argument("--peaks", help="Peaks BED/TSV (headerless). Must have >=3 columns: chr, start, end.")

    # Table mode (table -> table)
    p.add_argument("--table", help="Existing feature table TSV (must have chr/start/end).")

    p.add_argument("--phylop-bw", required=True, help="PhyloP bigWig in the TARGET assembly.")
    p.add_argument("--output", required=True, help="Output TSV path.")
    p.add_argument("--chain", help="UCSC chain file for source->target liftOver (optional).")
    p.add_argument("--liftover-exe", default="liftOver", help="Path to UCSC liftOver binary.")
    p.add_argument("--allow-missing-chroms", action="store_true", help="Missing chrom -> 0 instead of error.")
    p.add_argument("--max-len", type=int, default=5000, help="Intervals longer than this get 0.")
    p.add_argument("--out-col", default="phyloP_mean", help="Name of appended PhyloP column.")
    p.add_argument("--drop-lifted-coords", action="store_true", help="Do not keep lifted target coords.")
    p.add_argument("--quiet", action="store_true")
    return p


def main() -> None:
    args = build_parser().parse_args()

    if (args.peaks is None) == (args.table is None):
        raise SystemExit("Provide exactly one of --peaks OR --table")

    if args.peaks:
        peaks_df = load_peaks(args.peaks, quiet=args.quiet)

        # Optionally liftOver before scoring
        if args.chain:
            peaks_df = _run_liftover(
                peaks_df, chain_file=args.chain, liftover_exe=args.liftover_exe, quiet=args.quiet
            )
            chr_c, start_c, end_c = "chr_target", "start_target", "end_target"
        else:
            chr_c, start_c, end_c = "chr", "start", "end"

        log(f"Opening PhyloP bigWig: {args.phylop_bw}", args.quiet)
        bw = pyBigWig.open(args.phylop_bw)

        scores = []
        missing = 0
        too_long = 0
        for _, row in peaks_df.iterrows():
            chrom = row.get(chr_c)
            if pd.isna(chrom):
                missing += 1
                scores.append(0.0)
                continue

            start = int(row[start_c])
            end = int(row[end_c])

            if end <= start:
                scores.append(0.0)
                continue

            if (end - start) > args.max_len:
                too_long += 1
                scores.append(0.0)
                continue

            try:
                scores.append(_mean_bigwig(bw, str(chrom), start, end))
            except RuntimeError:
                missing += 1
                if args.allow_missing_chroms:
                    scores.append(0.0)
                else:
                    raise

        bw.close()

        if missing and not args.quiet:
            log(f"Warning: {missing} peaks missing target coordinates or chrom not in bigWig; filled as 0.", args.quiet)
        if too_long and not args.quiet:
            log(f"Warning: {too_long} peaks exceeded max_len={args.max_len}; filled as 0.", args.quiet)

        out_df = peaks_df.copy()
        out_df[args.out_col] = scores
        if args.chain and args.drop_lifted_coords:
            out_df = out_df.drop(columns=["chr_target", "start_target", "end_target"], errors="ignore")

        ensure_parent_dir(args.output)
        out_df.to_csv(args.output, sep="\t", index=False)
        log(f"Wrote PhyloP table: {args.output}", args.quiet)

    else:
        add_phylop(
            table_tsv=args.table,
            phylop_bw=args.phylop_bw,
            out_tsv=args.output,
            chain_file=args.chain,
            liftover_exe=args.liftover_exe,
            allow_missing_chroms=args.allow_missing_chroms,
            max_len=args.max_len,
            quiet=args.quiet,
            out_col=args.out_col,
            keep_lifted_coords=not args.drop_lifted_coords,
        )


if __name__ == "__main__":
    main()