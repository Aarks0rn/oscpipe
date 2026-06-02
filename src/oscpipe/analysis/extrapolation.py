"""Extrapolate an oligomer property to the polymer (infinite-chain) limit.

Conjugated-oligomer properties (optical gap, HOMO, …) vary approximately
linearly in 1/n, where n is the number of repeat units. A linear fit of
value vs 1/n gives the polymer limit as the intercept at 1/n -> 0. The fit
quality (r_squared) flags when more / longer oligomers are needed before the
limit can be trusted. Pure function — no IO, no DFT.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np


@dataclass
class ExtrapolationResult:
    limit: float  # value at 1/n -> 0 (polymer limit), i.e. the fit intercept
    slope: float  # d(value)/d(1/n)
    r_squared: float  # linear-fit quality in [0, 1]
    n_points: int


def extrapolate_inverse_n(points: Sequence[tuple[int, float]]) -> ExtrapolationResult:
    """Fit ``value`` vs ``1/n`` linearly and return the intercept (polymer limit).

    ``points`` is a sequence of ``(n, value)`` with n >= 1. Raises ``ValueError``
    for fewer than two points or any ``n < 1``.
    """
    if len(points) < 2:
        raise ValueError(f"need >= 2 oligomer points to extrapolate, got {len(points)}")
    if any(n < 1 for n, _ in points):
        raise ValueError("all oligomer lengths n must be >= 1")

    x = np.array([1.0 / n for n, _ in points])
    y = np.array([v for _, v in points])
    slope, intercept = np.polyfit(x, y, 1)

    residual = y - (slope * x + intercept)
    ss_res = float(np.sum(residual**2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0

    return ExtrapolationResult(
        limit=float(intercept),
        slope=float(slope),
        r_squared=float(r_squared),
        n_points=len(points),
    )
