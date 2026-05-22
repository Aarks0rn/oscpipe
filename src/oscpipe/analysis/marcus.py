"""Marcus rate + λ_reorg (4-point Nelsen scheme).

This is the single implementation. The OSC-pipeline archive had two
(`surface/marcus_ct_rate.py` and `computational/pipeline.py:compute_marcus_rate`)
— do not reintroduce a second one.

The 4-point workflow:
    1. neutral geometry at neutral charge   (E_n_n)
    2. neutral geometry at cation charge    (E_n_c)
    3. cation geometry at cation charge     (E_c_c)
    4. cation geometry at neutral charge    (E_c_n)
    λ_hole = (E_n_c - E_n_n) + (E_c_n - E_c_c)
"""

from __future__ import annotations

import math

_KB_EV = 8.617333262e-5  # eV/K
_HBAR_EVS = 6.582119569e-16  # eV·s


def lambda_hole_from_4_points(e_n_n: float, e_n_c: float, e_c_c: float, e_c_n: float) -> float:
    """All energies in eV, returns λ_hole in eV."""
    return (e_n_c - e_n_n) + (e_c_n - e_c_c)


def marcus_rate(
    lambda_ev: float, j_ev: float, delta_g_ev: float, temperature_k: float = 300.0
) -> float:
    """Classical Marcus rate constant in s^-1.

    k = (2π J² / ℏ) × (1 / √(4π λ k_B T)) × exp(−(ΔG + λ)² / (4λ k_B T))
    """
    kT = _KB_EV * temperature_k
    prefactor = (2.0 * math.pi * j_ev**2 / _HBAR_EVS) / math.sqrt(4.0 * math.pi * lambda_ev * kT)
    exponent = -((delta_g_ev + lambda_ev) ** 2) / (4.0 * lambda_ev * kT)
    return prefactor * math.exp(exponent)
