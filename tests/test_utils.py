# tests/test_utils.py
"""
Tests for pypeakranker._utils helper functions.
"""

import pandas as pd
import pytest

from pypeakranker._utils import load_peaks, log, ensure_parent_dir


def test_load_peaks_basic(tmp_path):
    """Basic 3-column BED file loads correctly."""
    bed = tmp_path / "peaks.bed"
    bed.write_text("chr1\t100\t200\nchr2\t300\t400\n")

    df = load_peaks(str(bed), quiet=True)
    assert list(df.columns[:3]) == ["chr", "start", "end"]
    assert len(df) == 2
    assert df["start"].dtype == int
    assert df["end"].dtype == int


def test_load_peaks_deduplicates(tmp_path):
    """Duplicate peaks are removed."""
    bed = tmp_path / "peaks.bed"
    bed.write_text("chr1\t100\t200\nchr1\t100\t200\nchr1\t300\t400\n")

    df = load_peaks(str(bed), quiet=True)
    assert len(df) == 2


def test_load_peaks_drops_inverted(tmp_path):
    """Peaks with end <= start are dropped."""
    bed = tmp_path / "peaks.bed"
    bed.write_text("chr1\t200\t100\nchr1\t300\t400\n")

    df = load_peaks(str(bed), quiet=True)
    assert len(df) == 1
    assert df.iloc[0]["start"] == 300


def test_load_peaks_drops_nonnumeric(tmp_path):
    """Rows with non-numeric start/end are dropped."""
    bed = tmp_path / "peaks.bed"
    bed.write_text("chr1\tabc\t200\nchr1\t300\t400\n")

    df = load_peaks(str(bed), quiet=True)
    assert len(df) == 1


def test_load_peaks_too_few_columns(tmp_path):
    """File with <3 columns raises ValueError."""
    bed = tmp_path / "peaks.bed"
    bed.write_text("chr1\t100\n")

    with pytest.raises(ValueError, match="at least 3 columns"):
        load_peaks(str(bed), quiet=True)


def test_load_peaks_extra_columns(tmp_path):
    """Extra columns beyond chr/start/end are preserved."""
    bed = tmp_path / "peaks.bed"
    bed.write_text("chr1\t100\t200\tpeak1\t500\n")

    df = load_peaks(str(bed), quiet=True)
    assert "col4" in df.columns
    assert "col5" in df.columns


def test_ensure_parent_dir(tmp_path):
    """Parent directories are created."""
    path = str(tmp_path / "a" / "b" / "output.tsv")
    ensure_parent_dir(path)
    import os
    assert os.path.isdir(str(tmp_path / "a" / "b"))


def test_log_quiet(capsys):
    """log() suppresses output when quiet=True."""
    log("hello", quiet=True)
    assert capsys.readouterr().out == ""


def test_log_verbose(capsys):
    """log() prints when quiet=False."""
    log("hello", quiet=False)
    assert "hello" in capsys.readouterr().out
