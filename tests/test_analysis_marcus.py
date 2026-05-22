"""Tests for oscpipe.analysis.marcus — lambda_hole_from_4_points and marcus_rate."""

import math

import pytest

from oscpipe.analysis.marcus import lambda_hole_from_4_points, marcus_rate


def test_lambda_hole_symmetric():
    # Both half-terms equal → λ = 2 × one term
    assert lambda_hole_from_4_points(0.0, 0.1, 0.0, 0.1) == pytest.approx(0.2)


def test_lambda_hole_asymmetric():
    # (0.05 - 0.0) + (0.08 - 0.03) = 0.05 + 0.05 = 0.10
    assert lambda_hole_from_4_points(0.0, 0.05, 0.03, 0.08) == pytest.approx(0.10)


def test_lambda_hole_zero():
    assert lambda_hole_from_4_points(1.0, 1.0, 1.0, 1.0) == pytest.approx(0.0)


def test_marcus_rate_known_point():
    # λ=0.3 eV, J=0.05 eV, ΔG=0, T=300 K → 4.2014e12 s⁻¹
    k = marcus_rate(lambda_ev=0.3, j_ev=0.05, delta_g_ev=0.0, temperature_k=300.0)
    assert k == pytest.approx(4.2014e12, rel=1e-3)


def test_marcus_rate_positive():
    k = marcus_rate(0.2, 0.03, 0.0, 300.0)
    assert k > 0


def test_marcus_rate_symmetric_maximum():
    # Rate is maximised at ΔG = -λ (activationless regime).
    lam = 0.3
    k_max = marcus_rate(lam, 0.05, delta_g_ev=-lam)
    k_off = marcus_rate(lam, 0.05, delta_g_ev=0.0)
    assert k_max > k_off


def test_marcus_rate_temperature_dependence():
    # Higher T → faster rate for normal (non-inverted) region.
    k_low = marcus_rate(0.3, 0.05, 0.0, temperature_k=200.0)
    k_high = marcus_rate(0.3, 0.05, 0.0, temperature_k=400.0)
    assert k_high > k_low


def test_marcus_rate_j_scaling():
    # Rate ∝ J² — doubling J quadruples rate.
    k1 = marcus_rate(0.3, 0.05, 0.0)
    k2 = marcus_rate(0.3, 0.10, 0.0)
    assert k2 == pytest.approx(4.0 * k1, rel=1e-6)


def test_marcus_rate_default_temperature():
    k_default = marcus_rate(0.3, 0.05, 0.0)
    k_300 = marcus_rate(0.3, 0.05, 0.0, temperature_k=300.0)
    assert k_default == pytest.approx(k_300)
