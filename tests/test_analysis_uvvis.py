"""Tests for `oscpipe.analysis.uvvis.broaden`."""

from __future__ import annotations

import numpy as np

from oscpipe.analysis.uvvis import _HC_NM_EV, broaden
from oscpipe.dft.gaussian import ExcitedState


def test_broaden_peak_position_matches_state_wavelength():
    e = 3.0  # eV
    wl_expected = _HC_NM_EV / e
    states = [ExcitedState(n=1, energy_ev=e, wavelength_nm=wl_expected, oscillator_strength=1.0)]
    wl, intensity = broaden(states, sigma_ev=0.1)
    peak_wl = wl[int(np.argmax(intensity))]
    assert abs(peak_wl - wl_expected) < 2.0  # within grid resolution


def test_broaden_accepts_dict_inputs():
    states = [{"n": 1, "energy_ev": 4.0, "wavelength_nm": 309.96, "f": 0.5}]
    wl, intensity = broaden(states)
    assert wl.shape == intensity.shape
    assert intensity.max() > 0


def test_broaden_two_states_sum_at_overlap():
    states = [
        ExcitedState(n=1, energy_ev=3.0, wavelength_nm=413.3, oscillator_strength=1.0),
        ExcitedState(n=2, energy_ev=3.0, wavelength_nm=413.3, oscillator_strength=2.0),
    ]
    _, intensity = broaden(states, sigma_ev=0.2)
    assert abs(intensity.max() - 3.0) < 1e-3


def test_broaden_zero_oscillator_gives_zero_spectrum():
    states = [ExcitedState(n=1, energy_ev=3.0, wavelength_nm=413.3, oscillator_strength=0.0)]
    _, intensity = broaden(states)
    assert intensity.max() == 0.0
