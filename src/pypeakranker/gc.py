#!/usr/bin/env python3
"""
gc.py (pypeakranker)

Compute GC fraction per peak from a reference genome FASTA.


- add_gc(table_tsv, reference_fasta, out_tsv): read an existing feature table (chr/start/end),
  append GC_content, write out (can overwrite same file).

Also supports standalone peaks -> gc table:
python -m pypeakranker.gc --peaks ... --reference-fasta ... --output ...
"""

from __future__ import annotations

import argparse
import os

import pandas as pd
from pyfaidx import Fasta
from pypeakranker._utils import log, ensure_parent_dir, load_peaks


def gc_fraction(seq: str) -> float:
    s = seq.upper()
    if not s:
        return float("nan")
    gc = s.count("G") + s.count("C")
    return gc / len(s)


# -------------------------
# Package API
# -------------------------

def add_gc(
    table_tsv: str,
    reference_fasta: str,
    out_tsv: str,
    allow_missing_chroms: bool = False,
    quiet: bool = False,
) -> None:
    """
    Read an existing feature table TSV (must contain chr/start/end),
    append GC_content, and write out.
    """
    base = pd.read_csv(table_tsv, sep="\t")
    needed = {"chr", "start", "end"}
    if not needed.issubset(base.columns):
        raise ValueError(f"Table must contain columns {needed}")

    if base.duplicated(["chr", "start", "end"]).any():
        raise ValueError("Table has duplicate (chr,start,end) rows; cannot merge safely.")

    log(f"Loading FASTA: {reference_fasta}", quiet)
    fasta = Fasta(reference_fasta)

    gc_vals = []
    missing = 0
    for _, row in base.iterrows():
        chrom = str(row["chr"])
        start = int(row["start"])
        end = int(row["end"])
        try:
            seq = fasta[chrom][start:end].seq
            gc_vals.append(gc_fraction(seq))
        except KeyError:
            missing += 1
            if allow_missing_chroms:
                gc_vals.append(pd.NA)
            else:
                raise KeyError(
                    f"Chromosome '{chrom}' not found in FASTA. "
                    f"Re-run with --allow-missing-chroms to fill as NA."
                )
        except Exception:
            gc_vals.append(pd.NA)

    if missing and not quiet:
        log(f"Warning: {missing} peaks had chromosomes not found in FASTA.", quiet)

    out = base.copy()
    out["GC_content"] = gc_vals

    ensure_parent_dir(out_tsv)
    out.to_csv(out_tsv, sep="\t", index=False)
    log(f"Wrote table with GC_content: {out_tsv}", quiet)


# -------------------------
# Standalone script CLI
# -------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Compute GC fraction per peak from a reference genome FASTA.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Backward-compatible standalone mode (peaks -> gc table)
    p.add_argument("--peaks", help="Peaks BED/TSV (headerless). Must have >=3 columns: chr, start, end.")

    # Option B mode (table -> table)
    p.add_argument("--table", help="Existing feature table TSV (must have chr/start/end).")

    p.add_argument("--reference-fasta", required=True, help="Reference genome FASTA.")
    p.add_argument("--output", required=True, help="Output TSV path.")
    p.add_argument("--allow-missing-chroms", action="store_true")
    p.add_argument("--quiet", action="store_true")
    return p


def main() -> None:
    args = build_parser().parse_args()

    if (args.peaks is None) == (args.table is None):
        raise SystemExit("Provide exactly one of --peaks OR --table")

    if args.peaks:
        peaks_df = load_peaks(args.peaks, quiet=args.quiet)
        log(f"Loading FASTA: {args.reference_fasta}", args.quiet)
        fasta = Fasta(args.reference_fasta)

        gc_vals = []
        missing = 0
        for _, row in peaks_df.iterrows():
            chrom = str(row["chr"])
            start = int(row["start"])
            end = int(row["end"])
            try:
                seq = fasta[chrom][start:end].seq
                gc_vals.append(gc_fraction(seq))
            except KeyError:
                missing += 1
                if args.allow_missing_chroms:
                    gc_vals.append(pd.NA)
                else:
                    raise
            except Exception:
                gc_vals.append(pd.NA)

        if missing and not args.quiet:
            log(f"Warning: {missing} peaks had chromosomes not found in FASTA.", args.quiet)

        out_df = peaks_df.copy()
        out_df["GC_content"] = gc_vals
        ensure_parent_dir(args.output)
        out_df.to_csv(args.output, sep="\t", index=False)
        log(f"Wrote GC table: {args.output}", args.quiet)

    else:
        add_gc(
            table_tsv=args.table,
            reference_fasta=args.reference_fasta,
            out_tsv=args.output,
            allow_missing_chroms=args.allow_missing_chroms,
            quiet=args.quiet,
        )


if __name__ == "__main__":
    main()