# tests/test_gc.py
"""
Tests for pypeakranker.gc module.
"""

import pandas as pd
import pytest

from pypeakranker.gc import gc_fraction, add_gc
from pypeakranker.summary import init_table


class TestGcFraction:
    def test_all_gc(self):
        assert gc_fraction("GCGCGC") == pytest.approx(1.0)

    def test_all_at(self):
        assert gc_fraction("ATATAT") == pytest.approx(0.0)

    def test_mixed(self):
        assert gc_fraction("ACGT") == pytest.approx(0.5)

    def test_lowercase(self):
        assert gc_fraction("acgt") == pytest.approx(0.5)

    def test_empty_string(self):
        import math

        assert math.isnan(gc_fraction(""))


class TestAddGc:
    def test_basic_pipeline(self, tmp_path):
        """init → add_gc produces correct GC values."""
        peaks = tmp_path / "peaks.bed"
        peaks.write_text("chr1\t0\t10\nchr1\t10\t20\n")

        fasta = tmp_path / "genome.fa"
        fasta.write_text(">chr1\nACGTACGTACGTACGTACGTACGTACGTACGT\n")

        table = tmp_path / "features.tsv"
        init_table(str(peaks), str(table), quiet=True)

        add_gc(
            table_tsv=str(table),
            reference_fasta=str(fasta),
            out_tsv=str(table),
            allow_missing_chroms=True,
            quiet=True,
        )

        df = pd.read_csv(table, sep="\t")
        assert "GC_content" in df.columns
        assert df["GC_content"].between(0, 1).all()

    def test_missing_chrom_raises(self, tmp_path):
        """Missing chrom without allow_missing_chroms raises KeyError."""
        peaks = tmp_path / "peaks.bed"
        peaks.write_text("chrX\t0\t10\n")

        fasta = tmp_path / "genome.fa"
        fasta.write_text(">chr1\nACGTACGTACGTACGTACGT\n")

        table = tmp_path / "features.tsv"
        init_table(str(peaks), str(table), quiet=True)

        with pytest.raises(KeyError):
            add_gc(
                table_tsv=str(table),
                reference_fasta=str(fasta),
                out_tsv=str(tmp_path / "out.tsv"),
                allow_missing_chroms=False,
                quiet=True,
            )

    def test_missing_chrom_allowed(self, tmp_path):
        """Missing chrom with allow_missing_chroms fills NA."""
        peaks = tmp_path / "peaks.bed"
        peaks.write_text("chrX\t0\t10\n")

        fasta = tmp_path / "genome.fa"
        fasta.write_text(">chr1\nACGTACGTACGTACGTACGT\n")

        table = tmp_path / "features.tsv"
        init_table(str(peaks), str(table), quiet=True)

        add_gc(
            table_tsv=str(table),
            reference_fasta=str(fasta),
            out_tsv=str(tmp_path / "out.tsv"),
            allow_missing_chroms=True,
            quiet=True,
        )

        df = pd.read_csv(tmp_path / "out.tsv", sep="\t")
        assert df["GC_content"].isna().all()
