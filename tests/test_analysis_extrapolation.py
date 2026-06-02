"""Tests for analysis/extrapolation: fit a property vs 1/n to the polymer limit.

Conjugated-oligomer properties (optical gap, HOMO) vary ~linearly in 1/n; the
intercept at 1/n -> 0 is the polymer (infinite-chain) limit. Pure function.
"""

import pytest

from oscpipe.analysis.extrapolation import extrapolate_inverse_n


def test_perfect_linear_recovers_intercept_and_slope():
    # value = 2.0 + 3.0 * (1/n)  =>  n=1:5.0, n=2:3.5, n=3:3.0
    res = extrapolate_inverse_n([(1, 5.0), (2, 3.5), (3, 3.0)])
    assert res.limit == pytest.approx(2.0)
    assert res.slope == pytest.approx(3.0)
    assert res.r_squared == pytest.approx(1.0)
    assert res.n_points == 3


def test_decreasing_gap_extrapolates_below_shortest_oligomer():
    # gap shrinks with length; polymer limit sits below the longest computed.
    res = extrapolate_inverse_n([(1, 3.0), (2, 2.4), (3, 2.2)])
    assert res.limit < 2.2
    assert res.slope > 0  # value rises with 1/n (shorter chain => larger gap)


def test_requires_at_least_two_points():
    with pytest.raises(ValueError):
        extrapolate_inverse_n([(1, 3.0)])


def test_rejects_n_below_one():
    with pytest.raises(ValueError):
        extrapolate_inverse_n([(0, 3.0), (2, 2.4)])
