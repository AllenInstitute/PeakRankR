# tests/test_smoke.py

import pandas as pd
import pytest

from pypeakranker.summary import init_table, add_signal
from pypeakranker.gc import add_gc


def test_pipeline_init_and_gc(tmp_path):
    peaks = tmp_path / "peaks.bed"
    peaks.write_text("chr1\t0\t10\nchr1\t10\t20\n")

    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nACGTACGTACGTACGTACGTACGTACGTACGT\n")

    table = tmp_path / "features.tsv"
    init_table(str(peaks), str(table))

    add_gc(
        table_tsv=str(table),
        reference_fasta=str(fasta),
        out_tsv=str(table),
        allow_missing_chroms=True,
    )

    df = pd.read_csv(table, sep="\t")

    assert {"chr", "start", "end", "GC_content"}.issubset(df.columns)
    assert len(df) == 2
    assert df["GC_content"].between(0, 1).all()
    assert df["GC_content"].notna().all()
    assert (df["GC_content"] == 0.5).all()


def test_public_api_imports():
    import pypeakranker

    assert hasattr(pypeakranker, "init_table")
    assert hasattr(pypeakranker, "add_signal")
    assert hasattr(pypeakranker, "add_gc")
    assert hasattr(pypeakranker, "add_phylop")
    assert hasattr(pypeakranker, "add_moments")
    assert hasattr(pypeakranker, "rank_by_specificity")


def test_add_signal_missing_bigwig_raises(tmp_path):
    table = tmp_path / "table.tsv"
    table.write_text("chr\tstart\tend\nchr1\t0\t10\n")

    # pyBigWig.open raises RuntimeError when it can't open a file
    with pytest.raises((FileNotFoundError, OSError, RuntimeError)):
        add_signal(
            table_tsv=str(table),
            bigwig_files=[str(tmp_path / "nonexistent.bw")],
            out_tsv=str(tmp_path / "out.tsv"),
        )

# ── rank_by_specificity ───────────────────────────────────────────────────

def test_rank_by_specificity_basic(tmp_path):
    """Rank peaks by specificity — peak with highest target/bg ratio gets rank 1."""
    table = tmp_path / "features.tsv"
    # Three peaks, two groups (target=A, background=A+B)
    # Peak 0: A=1.0, B=0.0 → ratio=1/0.5=2.0 (most specific to A)
    # Peak 1: A=0.5, B=0.5 → ratio=0.5/0.5=1.0
    # Peak 2: A=0.2, B=0.8 → ratio=0.2/0.5=0.4 (least specific to A)
    table.write_text(
        "chr\tstart\tend\tA_mean\tB_mean\n"
        "chr1\t0\t100\t1.0\t0.0\n"
        "chr1\t100\t200\t0.5\t0.5\n"
        "chr1\t200\t300\t0.2\t0.8\n"
    )

    from pypeakranker.rank import rank_by_specificity
    out = rank_by_specificity(
        table_tsv=str(table),
        target_cols=["A_mean"],
        background_cols=["A_mean", "B_mean"],
        out_tsv=str(tmp_path / "ranked.tsv"),
        quiet=True,
    )

    assert "specificity_score" in out.columns
    assert "specificity_rank" in out.columns
    assert out["specificity_score"].between(0, 1).all()
    # Rank 1 should be the peak most specific to A
    top = out[out["specificity_rank"] == 1].iloc[0]
    assert top["start"] == 0   # peak 0 is most specific to A
    assert top["specificity_score"] == pytest.approx(1.0)


def test_rank_by_specificity_uniform_signal(tmp_path):
    """When all peaks have equal signal, all specificity scores should be 0."""
    table = tmp_path / "features.tsv"
    table.write_text(
        "chr\tstart\tend\tA_mean\tB_mean\n"
        "chr1\t0\t100\t1.0\t1.0\n"
        "chr1\t100\t200\t1.0\t1.0\n"
    )
    from pypeakranker.rank import rank_by_specificity
    out = rank_by_specificity(
        table_tsv=str(table),
        target_cols=["A_mean"],
        background_cols=["A_mean", "B_mean"],
        out_tsv=str(tmp_path / "ranked.tsv"),
        quiet=True,
    )
    assert (out["specificity_score"] == 0).all()


def test_rank_by_specificity_missing_column_raises(tmp_path):
    """Asking for a column not in the table should raise ValueError."""
    table = tmp_path / "features.tsv"
    table.write_text(
        "chr\tstart\tend\tA_mean\n"
        "chr1\t0\t100\t1.0\n"
    )
    from pypeakranker.rank import rank_by_specificity
    with pytest.raises(ValueError, match="not found in table"):
        rank_by_specificity(
            table_tsv=str(table),
            target_cols=["NONEXISTENT"],
            background_cols=["A_mean"],
            out_tsv=str(tmp_path / "ranked.tsv"),
            quiet=True,
        )


def test_rank_by_specificity_zero_background(tmp_path):
    """Peaks where all groups have zero signal should get specificity_score=0."""
    table = tmp_path / "features.tsv"
    table.write_text(
        "chr\tstart\tend\tA_mean\tB_mean\n"
        "chr1\t0\t100\t0.0\t0.0\n"
        "chr1\t100\t200\t1.0\t0.0\n"
    )
    from pypeakranker.rank import rank_by_specificity
    out = rank_by_specificity(
        table_tsv=str(table),
        target_cols=["A_mean"],
        background_cols=["A_mean", "B_mean"],
        out_tsv=str(tmp_path / "ranked.tsv"),
        quiet=True,
    )
    zero_peak = out[out["start"] == 0]
    assert zero_peak["specificity_score"].values[0] == pytest.approx(0.0)


def test_rank_by_specificity_multiple_target_cols(tmp_path):
    """Multiple target columns should be averaged before computing specificity."""
    table = tmp_path / "features.tsv"
    table.write_text(
        "chr\tstart\tend\tA1_mean\tA2_mean\tB_mean\n"
        "chr1\t0\t100\t1.0\t0.8\t0.1\n"
        "chr1\t100\t200\t0.2\t0.1\t0.9\n"
    )
    from pypeakranker.rank import rank_by_specificity
    out = rank_by_specificity(
        table_tsv=str(table),
        target_cols=["A1_mean", "A2_mean"],
        background_cols=["A1_mean", "A2_mean", "B_mean"],
        out_tsv=str(tmp_path / "ranked.tsv"),
        quiet=True,
    )
    # Peak 0 should be more specific to A (target) than peak 1
    rank_peak0 = out[out["start"] == 0]["specificity_rank"].values[0]
    rank_peak1 = out[out["start"] == 100]["specificity_rank"].values[0]
    assert rank_peak0 < rank_peak1
