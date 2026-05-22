"""Tests for `oscpipe.analysis.indo.transfer_integral`."""

from __future__ import annotations

import math

from oscpipe.analysis.indo import transfer_integral


def _synth_dimer_log(occ_hartree: list[float], virt_hartree: list[float]) -> str:
    """Minimal Gaussian log with controllable alpha occ/virt eigenvalues."""
    occ_line = " Alpha  occ. eigenvalues -- " + "".join(f"{v:>10.5f}" for v in occ_hartree)
    virt_line = " Alpha virt. eigenvalues -- " + "".join(f"{v:>10.5f}" for v in virt_hartree)
    return (
        " SCF Done:  E(RB3LYP) =       -2.34000000000 A.U. after   8 cycles\n"
        + occ_line
        + "\n"
        + virt_line
        + "\n"
        " Normal termination of Gaussian 16 at ...\n"
    )


def test_transfer_integral_half_homo_splitting(tmp_path):
    # HOMO = -0.20 Ha, HOMO-1 = -0.21 Ha → ΔE = 0.01 Ha ≈ 0.272 eV → J = 0.136 eV
    log = tmp_path / "dimer.log"
    log.write_text(_synth_dimer_log(occ_hartree=[-0.30, -0.21, -0.20], virt_hartree=[0.05, 0.10]))
    j = transfer_integral(str(log))
    expected_ev = 0.01 * 27.211386245988 / 2.0
    assert math.isclose(j, expected_ev, rel_tol=1e-4)


def test_transfer_integral_returns_non_negative(tmp_path):
    log = tmp_path / "dimer.log"
    # Inverted order — magnitude still positive.
    log.write_text(_synth_dimer_log(occ_hartree=[-0.20, -0.21], virt_hartree=[0.05, 0.10]))
    j = transfer_integral(str(log))
    assert j >= 0.0
