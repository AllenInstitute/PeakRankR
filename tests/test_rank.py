# tests/test_rank.py

import pandas as pd
import pytest

from pypeakranker.rank import rank_by_specificity


def make_table(tmp_path, rows):
    """Write a minimal feature table TSV and return path."""
    df = pd.DataFrame(rows)
    p = tmp_path / "features.tsv"
    df.to_csv(p, sep="\t", index=False)
    return str(p)


def test_rank_basic(tmp_path):
    """Peak with highest target/background ratio should get rank 1."""
    rows = {
        "chr": ["chr1", "chr1", "chr1"],
        "start": [100, 500, 1000],
        "end": [200, 600, 1100],
        # Target group: A
        # Background: A, B, C
        "A_mean": [10.0, 1.0, 0.5],  # peak0 has strongest A signal
        "B_mean": [1.0, 1.0, 2.0],
        "C_mean": [1.0, 1.0, 2.0],
    }
    table = make_table(tmp_path, rows)
    out = str(tmp_path / "ranked.tsv")

    result = rank_by_specificity(
        table_tsv=table,
        target_cols=["A_mean"],
        background_cols=["A_mean", "B_mean", "C_mean"],
        out_tsv=out,
        quiet=True,
    )

    # peak0 has A=10 vs avg_bg=(10+1+1)/3=4 → ratio=2.5 — highest
    assert result.iloc[0]["specificity_rank"] == 1
    # Check first row is peak0 (chr1:100-200)
    assert result.iloc[0]["start"] == 100


def test_rank_score_range(tmp_path):
    """Specificity scores must be in [0, 1]."""
    rows = {
        "chr": ["chr1"] * 5,
        "start": [i * 100 for i in range(5)],
        "end": [i * 100 + 50 for i in range(5)],
        "A_mean": [5.0, 2.0, 8.0, 1.0, 4.0],
        "B_mean": [1.0, 3.0, 1.0, 1.0, 5.0],
    }
    table = make_table(tmp_path, rows)
    out = str(tmp_path / "ranked.tsv")

    result = rank_by_specificity(
        table_tsv=table,
        target_cols=["A_mean"],
        background_cols=["A_mean", "B_mean"],
        out_tsv=out,
        quiet=True,
    )

    scores = result["specificity_score"]
    assert (scores >= 0).all()
    assert (scores <= 1).all()
    assert scores.max() == pytest.approx(1.0)
    assert scores.min() == pytest.approx(0.0)


def test_rank_all_zero_signal(tmp_path):
    """All-zero signal should not raise; all scores = 0, all ranks = 1."""
    rows = {
        "chr": ["chr1", "chr1"],
        "start": [100, 200],
        "end": [200, 300],
        "A_mean": [0.0, 0.0],
        "B_mean": [0.0, 0.0],
    }
    table = make_table(tmp_path, rows)
    out = str(tmp_path / "ranked.tsv")

    result = rank_by_specificity(
        table_tsv=table,
        target_cols=["A_mean"],
        background_cols=["A_mean", "B_mean"],
        out_tsv=out,
        quiet=True,
    )

    assert (result["specificity_score"] == 0.0).all()
    # all tied → all rank 1
    assert (result["specificity_rank"] == 1).all()


def test_rank_missing_target_col_raises(tmp_path):
    """Non-existent target column should raise ValueError."""
    rows = {
        "chr": ["chr1"],
        "start": [100],
        "end": [200],
        "A_mean": [5.0],
    }
    table = make_table(tmp_path, rows)
    out = str(tmp_path / "out.tsv")

    with pytest.raises(ValueError, match="not found in table"):
        rank_by_specificity(
            table_tsv=table,
            target_cols=["NONEXISTENT"],
            background_cols=["A_mean"],
            out_tsv=out,
            quiet=True,
        )


def test_rank_multiple_target_cols_mean(tmp_path):
    """Multiple target columns aggregated by mean."""
    rows = {
        "chr": ["chr1", "chr1"],
        "start": [100, 500],
        "end": [200, 600],
        "A1_mean": [8.0, 1.0],
        "A2_mean": [6.0, 1.0],  # mean target for peak0 = 7.0
        "B_mean": [1.0, 5.0],
    }
    table = make_table(tmp_path, rows)
    out = str(tmp_path / "ranked.tsv")

    result = rank_by_specificity(
        table_tsv=table,
        target_cols=["A1_mean", "A2_mean"],
        background_cols=["A1_mean", "A2_mean", "B_mean"],
        out_tsv=out,
        target_agg="mean",
        quiet=True,
    )

    # peak0: mean_target=(8+6)/2=7, avg_bg=(8+6+1)/3=5 → ratio=1.4
    # peak1: mean_target=(1+1)/2=1, avg_bg=(1+1+5)/3=2.33 → ratio=0.43
    # peak0 should be rank 1
    assert result.iloc[0]["start"] == 100


def test_rank_output_file_written(tmp_path):
    """Output TSV should be written to disk."""
    rows = {
        "chr": ["chr1"],
        "start": [100],
        "end": [200],
        "A_mean": [5.0],
        "B_mean": [1.0],
    }
    table = make_table(tmp_path, rows)
    out = str(tmp_path / "subdir" / "ranked.tsv")

    rank_by_specificity(
        table_tsv=table,
        target_cols=["A_mean"],
        background_cols=["A_mean", "B_mean"],
        out_tsv=out,
        quiet=True,
    )

    import os

    assert os.path.exists(out)
    df = pd.read_csv(out, sep="\t")
    assert "specificity_score" in df.columns
    assert "specificity_rank" in df.columns


def test_rank_public_api():
    """rank_by_specificity is accessible from the package root."""
    import pypeakranker

    assert hasattr(pypeakranker, "rank_by_specificity")
