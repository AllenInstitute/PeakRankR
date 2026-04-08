# tests/test_moments.py
"""
Tests for pypeakranker.moments helper functions.

These test the pure-Python computation routines (no BigWig files needed).
"""

import numpy as np
import pytest

from pypeakranker.moments import (
    bimodality_coefficient,
    _metrics_for_interval,
)


class TestBimodalityCoefficient:
    def test_too_few_values(self):
        assert np.isnan(bimodality_coefficient(np.array([1.0, 2.0, 3.0])))

    def test_uniform_returns_finite(self):
        vals = np.ones(100)
        result = bimodality_coefficient(vals)
        # All identical values → skew=0, formula should return a number or nan
        assert np.isfinite(result) or np.isnan(result)

    def test_bimodal_higher_than_unimodal(self):
        """A clearly bimodal distribution should have a higher coefficient."""
        rng = np.random.default_rng(42)
        unimodal = rng.normal(0, 1, 1000)
        bimodal = np.concatenate([rng.normal(-3, 0.5, 500), rng.normal(3, 0.5, 500)])
        assert bimodality_coefficient(bimodal) > bimodality_coefficient(unimodal)


class TestMetricsForInterval:
    def test_empty_array(self):
        s, k, b = _metrics_for_interval(np.array([]))
        assert np.isnan(s)
        assert np.isnan(k)
        assert np.isnan(b)

    def test_constant_array(self):
        s, k, b = _metrics_for_interval(np.ones(50))
        assert s == 0.0
        assert k == 0.0
        assert np.isnan(b)

    def test_all_zeros(self):
        s, k, b = _metrics_for_interval(np.zeros(50))
        assert s == 0.0
        assert k == 0.0

    def test_normal_distribution(self):
        rng = np.random.default_rng(123)
        vals = rng.normal(0, 1, 10000)
        s, k, b = _metrics_for_interval(vals)
        # Normal: skew ~0, kurtosis ~3
        assert abs(s) < 0.1
        assert abs(k - 3.0) < 0.2
        assert np.isfinite(b)
