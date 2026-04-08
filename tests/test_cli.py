# tests/test_cli.py
"""
Tests for the CLI parser and dispatch.
These tests exercise argument parsing and the main() routing logic
without requiring real BigWig files (except where noted).
"""

import subprocess
import sys


def test_cli_help():
    """pypeakranker --help exits 0."""
    result = subprocess.run(
        [sys.executable, "-m", "pypeakranker.cli", "--help"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert "pypeakranker" in result.stdout


def test_cli_init_subcommand(tmp_path):
    """pypeakranker init creates a valid feature table."""
    peaks = tmp_path / "peaks.bed"
    peaks.write_text("chr1\t100\t200\nchr1\t300\t400\n")
    out = tmp_path / "features.tsv"

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "pypeakranker.cli",
            "init",
            "--peaks",
            str(peaks),
            "--out",
            str(out),
            "--quiet",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert out.exists()

    import pandas as pd

    df = pd.read_csv(out, sep="\t")
    assert list(df.columns[:3]) == ["chr", "start", "end"]
    assert len(df) == 2


def test_cli_rank_specificity_subcommand(tmp_path):
    """pypeakranker rank-specificity runs end-to-end via CLI."""
    table = tmp_path / "features.tsv"
    table.write_text(
        "chr\tstart\tend\tA_mean\tB_mean\n"
        "chr1\t0\t100\t10.0\t1.0\n"
        "chr1\t100\t200\t1.0\t10.0\n"
    )
    out = tmp_path / "ranked.tsv"

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "pypeakranker.cli",
            "rank-specificity",
            "--table",
            str(table),
            "--target-cols",
            "A_mean",
            "--background-cols",
            "A_mean",
            "B_mean",
            "--out",
            str(out),
            "--quiet",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert out.exists()

    import pandas as pd

    df = pd.read_csv(out, sep="\t")
    assert "specificity_score" in df.columns
    assert "specificity_rank" in df.columns
    # Peak with A=10, B=1 should be rank 1
    assert df.iloc[0]["specificity_rank"] == 1


def test_cli_missing_subcommand():
    """Running without a subcommand should exit non-zero."""
    result = subprocess.run(
        [sys.executable, "-m", "pypeakranker.cli"],
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
