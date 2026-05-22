"""dft/gaussian tests.

Pure-function checks of the .com writer and is_log_complete tail-read.
Parser tests use real fixtures from user@203.0.113.10 (HF/STO-3G for H2,
B3LYP/6-31G* for benzene) committed under tests/fixtures/.
"""

from pathlib import Path

import ase
import pytest

from oscpipe.dft import gaussian

FIXTURES = Path(__file__).parent / "fixtures"

# ── write_com_properties ───────────────────────────────────────────────────


def _h2():
    return ase.Atoms("H2", positions=[(0.0, 0.0, 0.0), (0.74, 0.0, 0.0)])


def test_write_com_properties_has_route_line():
    com = gaussian.write_com_properties(
        _h2(),
        method="b3lyp",
        basis="6-31g*",
        charge=0,
        mult=1,
        nproc=4,
        mem="4GB",
        label="h2_test",
        chk="h2_test.chk",
    )
    assert "%nprocshared=4" in com
    assert "%mem=4GB" in com
    assert "%chk=h2_test.chk" in com
    assert "#p b3lyp/6-31g* opt pop=full" in com.lower()


def test_write_com_properties_has_charge_mult_and_coords():
    com = gaussian.write_com_properties(_h2(), "b3lyp", "6-31g*", 0, 1, 4, "4GB", "h2", "h2.chk")
    lines = com.splitlines()
    # charge/mult line should appear; then two H lines with x y z.
    charge_idx = lines.index("0 1")
    assert lines[charge_idx + 1].split()[0] == "H"
    assert lines[charge_idx + 2].split()[0] == "H"
    # coordinates within tolerance
    second_h = lines[charge_idx + 2].split()
    assert abs(float(second_h[1]) - 0.74) < 1e-6


def test_write_com_properties_trailing_blank_line():
    com = gaussian.write_com_properties(_h2(), "b3lyp", "6-31g*", 0, 1, 4, "4GB", "h2", "h2.chk")
    # Gaussian requires the molecule spec to end with a blank line.
    assert com.endswith("\n\n")


# ── is_log_complete ────────────────────────────────────────────────────────


def test_is_log_complete_normal_termination(tmp_path: Path):
    log = tmp_path / "ok.log"
    log.write_text("some output\n Normal termination of Gaussian 16 at ...\n")
    assert gaussian.is_log_complete(str(log)) is True


def test_is_log_complete_missing_file(tmp_path: Path):
    assert gaussian.is_log_complete(str(tmp_path / "nope.log")) is False


def test_is_log_complete_no_marker(tmp_path: Path):
    log = tmp_path / "wip.log"
    log.write_text("Entering Link 1\n... still running\n")
    assert gaussian.is_log_complete(str(log)) is False


# ── parse_properties (real fixture logs from user@203.0.113.10) ───────────


def test_parse_properties_h2():
    """H2 / HF/STO-3G — single occupied MO; dipole = 0 by symmetry."""
    r = gaussian.parse_properties(str(FIXTURES / "h2.log"))
    assert abs(r.homo_ev - (-15.732)) < 1e-2
    assert abs(r.lumo_ev - 18.235) < 1e-2
    assert r.gap_ev > 0
    assert abs(r.gap_ev - (r.lumo_ev - r.homo_ev)) < 1e-6
    assert abs(r.dipole_debye - 0.0) < 1e-2
    assert r.energy_ev < 0


def test_parse_properties_benzene():
    """Benzene / B3LYP/6-31G* — HOMO~−6.7 eV, dipole = 0 (D6h)."""
    r = gaussian.parse_properties(str(FIXTURES / "benzene.log"))
    assert abs(r.homo_ev - (-6.693)) < 1e-2
    assert abs(r.lumo_ev - 0.077) < 1e-2
    assert abs(r.gap_ev - 6.770) < 1e-2
    assert abs(r.dipole_debye - 0.0) < 1e-2
    assert r.energy_ev < -6000


# ── parse_dimer_orbitals ───────────────────────────────────────────────────


def test_parse_dimer_orbitals_benzene_returns_four_values():
    """Benzene has degenerate HOMO/HOMO-1 and LUMO/LUMO+1 (D6h symmetry)."""
    homo, hm1, lumo, lp1 = gaussian.parse_dimer_orbitals(str(FIXTURES / "benzene.log"))
    assert abs(homo - (-6.693)) < 1e-2
    assert abs(hm1 - (-6.693)) < 1e-2  # degenerate pair
    assert abs(lumo - 0.077) < 1e-2
    assert abs(lp1 - 0.077) < 1e-2  # degenerate pair


def test_parse_dimer_orbitals_non_degenerate(tmp_path: Path):
    """Synthetic log with split HOMO/HOMO-1 — J_hole = |Δ|/2."""
    log = tmp_path / "dimer.log"
    # Two occ lines → 10 eigenvalues each (in Hartree); HOMO = -0.3, HOMO-1 = -0.32
    log.write_text(
        " Alpha  occ. eigenvalues -- "
        " -1.00000 -0.90000 -0.80000 -0.70000 -0.60000\n"
        " Alpha  occ. eigenvalues -- "
        " -0.50000 -0.40000 -0.35000 -0.32000 -0.30000\n"
        " Alpha virt. eigenvalues --  "
        "  0.10000  0.12000  0.20000  0.30000  0.40000\n"
        " Normal termination of Gaussian 16 at ...\n"
    )
    HARTREE_TO_EV = 27.2114
    homo, hm1, lumo, lp1 = gaussian.parse_dimer_orbitals(str(log))
    assert abs(homo - (-0.30000 * HARTREE_TO_EV)) < 1e-4
    assert abs(hm1 - (-0.32000 * HARTREE_TO_EV)) < 1e-4
    assert abs(lumo - (0.10000 * HARTREE_TO_EV)) < 1e-4
    assert abs(lp1 - (0.12000 * HARTREE_TO_EV)) < 1e-4
    j_hole = abs(homo - hm1) / 2
    assert abs(j_hole - (0.01 * HARTREE_TO_EV)) < 1e-4


def test_parse_dimer_orbitals_insufficient_orbitals_raises(tmp_path: Path):
    log = tmp_path / "monomer.log"
    log.write_text(
        " Alpha  occ. eigenvalues --  -0.50000\n"
        " Alpha virt. eigenvalues --   0.10000\n"
        " Normal termination of Gaussian 16 at ...\n"
    )
    with pytest.raises(ValueError, match="≥2 occ"):
        gaussian.parse_dimer_orbitals(str(log))
