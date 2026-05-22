"""UV-Vis spectrum from TDDFT excited states.

Reads a list of ExcitedState (energy_ev, oscillator_strength), Gaussian-broadens
in energy space, returns (wavelength_nm, intensity) arrays for plotting.
"""

from __future__ import annotations

from typing import Iterable

import numpy as np

_HC_NM_EV = 1239.841984  # h*c in nm·eV


def _state_to_pair(s) -> tuple[float, float]:
    """Accept ExcitedState dataclass or dict with energy_ev + oscillator_strength|f."""
    if hasattr(s, "energy_ev"):
        return float(s.energy_ev), float(s.oscillator_strength)
    e = float(s["energy_ev"])
    f = float(s.get("oscillator_strength", s.get("f")))
    return e, f


def broaden(
    states: Iterable,
    sigma_ev: float = 0.2,
    wavelength_min_nm: float = 200.0,
    wavelength_max_nm: float = 800.0,
    n_points: int = 601,
) -> tuple[np.ndarray, np.ndarray]:
    """Gaussian-broaden TDDFT sticks in energy space; sample on wavelength grid.

    intensity(E) = Σ f_i · exp(-(E - E_i)² / (2 σ²))
    """
    wl = np.linspace(wavelength_min_nm, wavelength_max_nm, n_points)
    energy = _HC_NM_EV / wl
    intensity = np.zeros_like(wl)
    for s in states:
        e_i, f_i = _state_to_pair(s)
        intensity += f_i * np.exp(-0.5 * ((energy - e_i) / sigma_ev) ** 2)
    return wl, intensity
