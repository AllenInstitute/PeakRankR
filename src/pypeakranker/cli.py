import argparse

from pypeakranker.summary import init_table, add_signal
from pypeakranker.gc import add_gc
from pypeakranker.phylop import add_phylop
from pypeakranker.moments import add_moments
from pypeakranker.rank import rank_by_specificity


def build_parser():
    p = argparse.ArgumentParser(prog="pypeakranker")
    sub = p.add_subparsers(dest="cmd", required=True)

    # -------------------------
    # init
    # -------------------------
    p_init = sub.add_parser("init", help="Initialize a feature table from peaks")
    p_init.add_argument(
        "--peaks", required=True, help="Peaks BED/TSV (>=3 cols, headerless)"
    )
    p_init.add_argument("--out", required=True, help="Output TSV (feature table)")
    p_init.add_argument("--quiet", action="store_true")

    # -------------------------
    # add-signal
    # -------------------------
    p_sig = sub.add_parser(
        "add-signal", help="Add BigWig signal summaries to feature table"
    )
    p_sig.add_argument("--table", required=True, help="Existing feature table TSV")
    p_sig.add_argument("--bigwig-files", nargs="+", required=True, help="BigWig files")
    p_sig.add_argument(
        "--out", required=True, help="Output TSV (can overwrite same file)"
    )
    p_sig.add_argument(
        "--stat",
        default="sum",
        choices=["sum", "mean", "max"],
        help="Summary statistic",
    )
    p_sig.add_argument(
        "--suffix", default="summary", help="Suffix for new columns (e.g., summary)"
    )
    p_sig.add_argument(
        "--sample-name-mode", default="stem", choices=["stem", "filename"]
    )
    p_sig.add_argument("--keep-nans", action="store_true")
    p_sig.add_argument(
        "--allow-missing-chroms",
        action="store_true",
        help="Missing chrom -> 0 instead of error",
    )
    p_sig.add_argument("--quiet", action="store_true")

    # -------------------------
    # add-gc
    # -------------------------
    p_gc = sub.add_parser("add-gc", help="Add GC content to feature table")
    p_gc.add_argument("--table", required=True, help="Existing feature table TSV")
    p_gc.add_argument("--reference-fasta", required=True, help="Reference FASTA")
    p_gc.add_argument(
        "--out", required=True, help="Output TSV (can overwrite same file)"
    )
    p_gc.add_argument(
        "--allow-missing-chroms",
        action="store_true",
        help="Missing chrom -> NA instead of error",
    )
    p_gc.add_argument("--quiet", action="store_true")

    # -------------------------
    # add-phylop
    # -------------------------
    p_phy = sub.add_parser("add-phylop", help="Add mean PhyloP score to feature table")
    p_phy.add_argument("--table", required=True, help="Existing feature table TSV")
    p_phy.add_argument(
        "--phylop-bw", required=True, help="PhyloP bigWig (TARGET assembly)"
    )
    p_phy.add_argument(
        "--out", required=True, help="Output TSV (can overwrite same file)"
    )
    p_phy.add_argument(
        "--chain", help="UCSC chain file for source->target liftOver (optional)"
    )
    p_phy.add_argument(
        "--liftover-exe", default="liftOver", help="Path to UCSC liftOver binary"
    )
    p_phy.add_argument(
        "--allow-missing-chroms",
        action="store_true",
        help="Missing chrom -> 0 instead of error",
    )
    p_phy.add_argument(
        "--max-len", type=int, default=5000, help="Intervals longer than this get 0"
    )
    p_phy.add_argument(
        "--out-col", default="phyloP_mean", help="Name of appended PhyloP column"
    )
    p_phy.add_argument(
        "--drop-lifted-coords",
        action="store_true",
        help="Drop chr_target/start_target/end_target",
    )
    p_phy.add_argument("--quiet", action="store_true")

    # -------------------------
    # add-moments
    # -------------------------
    p_mom = sub.add_parser(
        "add-moments",
        help="Add skewness/kurtosis/bimodality per BigWig to feature table",
    )
    p_mom.add_argument("--table", required=True, help="Existing feature table TSV")
    p_mom.add_argument("--bigwig-files", nargs="+", required=True, help="BigWig files")
    p_mom.add_argument(
        "--out", required=True, help="Output TSV (can overwrite same file)"
    )
    p_mom.add_argument(
        "--prefix", default="", help="Optional prefix for appended columns"
    )
    p_mom.add_argument(
        "--allow-missing-chroms",
        action="store_true",
        help="Missing chrom -> NaN metrics instead of error",
    )
    p_mom.add_argument("--quiet", action="store_true")

    # -------------------------
    # rank-specificity
    # -------------------------
    p_rank = sub.add_parser("rank-specificity", help="Rank peaks by ATAC specificity")
    p_rank.add_argument(
        "--table",
        required=True,
        help="Feature table TSV (must have chr/start/end and signal columns)",
    )
    p_rank.add_argument(
        "--target-cols",
        nargs="+",
        required=True,
        help="Signal column(s) for the target group",
    )
    p_rank.add_argument(
        "--background-cols",
        nargs="+",
        required=True,
        help="Signal columns for ALL groups including target",
    )
    p_rank.add_argument("--out", required=True, help="Output TSV path")
    p_rank.add_argument(
        "--target-agg",
        default="mean",
        choices=["mean", "sum"],
        help="How to aggregate multiple --target-cols (default: mean)",
    )
    p_rank.add_argument(
        "--specificity-col",
        default="specificity_score",
        help="Name for the specificity score column",
    )
    p_rank.add_argument(
        "--rank-col", default="specificity_rank", help="Name for the rank column"
    )
    p_rank.add_argument("--quiet", action="store_true")

    return p


def main():
    p = build_parser()
    args = p.parse_args()

    if args.cmd == "init":
        init_table(args.peaks, args.out, quiet=args.quiet)

    elif args.cmd == "add-signal":
        add_signal(
            table_tsv=args.table,
            bigwig_files=args.bigwig_files,
            out_tsv=args.out,
            stat=args.stat,
            suffix=args.suffix,
            sample_name_mode=args.sample_name_mode,
            keep_nans=args.keep_nans,
            allow_missing_chroms=args.allow_missing_chroms,
            quiet=args.quiet,
        )

    elif args.cmd == "add-gc":
        add_gc(
            table_tsv=args.table,
            reference_fasta=args.reference_fasta,
            out_tsv=args.out,
            allow_missing_chroms=args.allow_missing_chroms,
            quiet=args.quiet,
        )

    elif args.cmd == "add-phylop":
        add_phylop(
            table_tsv=args.table,
            phylop_bw=args.phylop_bw,
            out_tsv=args.out,
            chain_file=args.chain,
            liftover_exe=args.liftover_exe,
            allow_missing_chroms=args.allow_missing_chroms,
            max_len=args.max_len,
            quiet=args.quiet,
            out_col=args.out_col,
            keep_lifted_coords=not args.drop_lifted_coords,
        )

    elif args.cmd == "rank-specificity":
        rank_by_specificity(
            table_tsv=args.table,
            target_cols=args.target_cols,
            background_cols=args.background_cols,
            out_tsv=args.out,
            target_agg=args.target_agg,
            specificity_col=args.specificity_col,
            rank_col=args.rank_col,
            quiet=args.quiet,
        )

    elif args.cmd == "add-moments":
        add_moments(
            table_tsv=args.table,
            bigwig_files=args.bigwig_files,
            out_tsv=args.out,
            allow_missing_chroms=args.allow_missing_chroms,
            quiet=args.quiet,
            prefix=args.prefix,
        )


if __name__ == "__main__":
    main()
