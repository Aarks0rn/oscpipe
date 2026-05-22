"""Transfer integral via ZINDO dimer HOMO/HOMO-1 splitting.

Reuses dimer construction from `oscpipe.chem.geometry`. Cofacial Δz default
3.5 Å (matches OSC-pipeline validation benchmarks).

Energy-splitting method (Newton):
    J_hole = |E_HOMO − E_HOMO-1| / 2
"""

from __future__ import annotations

from oscpipe.dft.gaussian import parse_dimer_orbitals


def transfer_integral(dimer_log_path: str) -> float:
    """Return J_hole in eV from a ZINDO/DFT dimer log."""
    homo, homo_m1, _lumo, _lumo_p1 = parse_dimer_orbitals(dimer_log_path)
    return abs(homo - homo_m1) / 2.0
