"""Gaussian .com writer and .log parser. Pure functions, no IO orchestration.

Writers take a config and an ASE Atoms object, return string contents.
Parsers take a path to a completed .log, return a dataclass of results.
Dispatch + lifecycle live in `oscpipe.dispatch`; this module never calls
subprocess or paramiko.
"""

from __future__ import annotations

import os
from dataclasses import dataclass

HARTREE_TO_EV = 27.2114


@dataclass
class PropertiesResult:
    homo_ev: float
    lumo_ev: float
    gap_ev: float
    dipole_debye: float
    energy_ev: float


@dataclass
class ExcitedState:
    n: int
    energy_ev: float
    wavelength_nm: float
    oscillator_strength: float


# ── writers ────────────────────────────────────────────────────────────────


def write_com_properties(
    atoms,
    method: str,
    basis: str,
    charge: int,
    mult: int,
    nproc: int,
    mem: str,
    label: str,
    chk: str,
    opt_route: str = "opt",
) -> str:
    """Return the .com text for an opt + pop=full properties job.

    ``opt_route`` is the bare opt keyword. Default ``"opt"`` keeps the tight
    default convergence that λ_h needs (reorganisation energy is sensitive to the
    neutral/cation geometry difference, so those opts must NOT be loosened). The
    oligomer sweep overrides it with damping + loose convergence: floppy backbones
    (TVT thienylene-vinylene arms) give a near-flat PES where the default trust
    radius overshoots — FBT n=1 thrashed with Max Displacement ~0.7 (target 0.0018)
    and ran to the 462-step cap without converging (job 71248). HOMO/KS-gap/optical
    gap are insensitive to the last bit of geometry convergence, so ``loose`` (Max
    Force 0.0025) is scientifically fine there and ``maxstep=10`` damps the
    overshoot. (Internal coords kept; if a long-n opt later hits 'FormBX'/'Bend
    failed' from a near-linear angle, add ``cartesian`` to its ``opt_route`` too.)
    """
    return _format_com(
        atoms=atoms,
        method=method,
        basis=basis,
        charge=charge,
        mult=mult,
        nproc=nproc,
        mem=mem,
        label=label,
        chk=chk,
        # scf=(xqc,vshift=100): long conjugated oligomers have near-degenerate
        # frontier orbitals (Gap≈0.02-0.06 au) where plain DIIS oscillates.
        # vshift level-shifts the virtuals so DIIS converges cheaply (xqc alone
        # converged but fell back to quadratic-convergence EVERY opt step ≈ 5h/step,
        # impractical over an opt). xqc stays as the safety net if DIIS still stalls.
        # vshift=300 over-shifted (≫ gap): on some conformers DIIS stalled at
        # RMSDP≈1e-6, never reached the tight criterion, fell to xqc and hung >15h
        # on the FIRST SCF point (job f665f0c6, 0 SCF Done). vshift=100 (~2x the
        # gap) still aids DIIS but converges tight. Do NOT raise back to 300
        # (stalls) or drop to 0 (5h/step). No-op cost when DIIS already converges.
        keywords=f"{opt_route} pop=full scf=(xqc,vshift=100)",
    )


def write_com_preopt(
    atoms,
    charge: int,
    mult: int,
    nproc: int,
    mem: str,
    label: str,
    chk: str,
    method: str = "pm6",
) -> str:
    """Return the .com text for a cheap semi-empirical geometry pre-optimisation.

    Large conjugated oligomers (n=3) start far from their minimum, so a direct
    B3LYP/6-31G** opt thrashes — huge RMS displacement, ~14 h/step, effectively
    never converging (job 71235, 5 steps in 3 days). A PM6 pre-opt lands a
    near-physical geometry far more cheaply; feeding that into the B3LYP opt cuts
    it to a handful of steps. PM6 is semi-empirical: no basis set.

    Two robustness flags, both forced by the n=3 backbone:

    * ``scf=(xqc,vshift=300)`` — the near-degenerate frontier (Gap≈0.1 eV at n=3)
      is a property of the conjugated system, NOT the functional, so PM6's SCF
      hits the SAME DIIS oscillation B3LYP does; plain PM6 dies with 'Convergence
      failure' at l502 (job 71238).
    * ``opt=cartesian`` — as the long oligomer relaxes it passes through a
      near-linear (≈180°) bond angle, which makes Gaussian's redundant *internal*
      coordinates singular ('Bend failed' → 'FormBX had a problem', l103 error
      termination, job 71239 — after 34 good steps). Optimising in Cartesian
      coordinates sidesteps the B-matrix entirely. Cheap here: PM6 is ~15 s/step.
    * ``loose,maxstep=10,maxcycles=300`` — a pre-opt only needs a near-physical
      geometry to seed the B3LYP opt; it does NOT need a tight stationary point.
      Floppy TVT backbones (soft thienylene-vinylene + 2D side-arm torsions) give a
      near-flat PES: the default trust radius overshoots, so Max Force oscillates
      and never settles even at the ``loose`` 0.0025 target — PBDT-TVT-FBT n=1 hit
      'Number of steps exceeded, NStep=241' first at the tight default (job 71241),
      then again with ``loose`` alone, oscillating 0.003–0.007 just above 0.0025
      (job 71244). ``loose`` relaxes convergence ~10x; ``maxstep=10`` (0.1 Bohr cap)
      damps the overshoot so it settles; ``maxcycles=300`` gives the smaller steps
      room. Harmless to the result: the B3LYP opt tightens the geometry afterward.
    """
    return _format_com(
        atoms=atoms,
        method=method,
        basis="",  # semi-empirical: no basis set
        charge=charge,
        mult=mult,
        nproc=nproc,
        mem=mem,
        label=label,
        chk=chk,
        keywords="opt=(cartesian,loose,maxstep=10,maxcycles=300) scf=(xqc,vshift=300)",  # see docstring
    )


def write_com_sp(
    atoms,
    method: str,
    basis: str,
    charge: int,
    mult: int,
    nproc: int,
    mem: str,
    label: str,
    chk: str,
) -> str:
    """Return the .com text for a single-point properties job (no opt).

    Used by the λ_reorg workflow to evaluate energies at a fixed geometry
    (the optimised neutral or cation structure from a sibling job).
    """
    return _format_com(
        atoms=atoms,
        method=method,
        basis=basis,
        charge=charge,
        mult=mult,
        nproc=nproc,
        mem=mem,
        label=label,
        chk=chk,
        keywords="pop=full scf=(xqc,vshift=300)",  # robust SCF — see write_com_properties
    )


def write_com_tddft(
    atoms,
    method: str,
    basis: str,
    charge: int,
    mult: int,
    nstates: int,
    nproc: int,
    mem: str,
    label: str,
    chk: str,
) -> str:
    """Return the .com text for a TDDFT excited-states job."""
    return _format_com(
        atoms=atoms,
        method=method,
        basis=basis,
        charge=charge,
        mult=mult,
        nproc=nproc,
        mem=mem,
        label=label,
        chk=chk,
        keywords=f"td=(nstates={nstates}) scf=(xqc,vshift=100)",  # robust SCF — see write_com_properties
    )


def _format_com(
    *,
    atoms,
    method: str,
    basis: str,
    charge: int,
    mult: int,
    nproc: int,
    mem: str,
    label: str,
    chk: str,
    keywords: str,
) -> str:
    # Semi-empirical methods (e.g. PM6) take no basis set — emit a bare method.
    route = f"{method}/{basis}" if basis else method
    header = [
        f"%nprocshared={nproc}",
        f"%mem={mem}",
        f"%chk={chk}",
        f"#p {route} {keywords}",
        "",
        label,
        "",
        f"{charge} {mult}",
    ]
    body = [
        f"{sym:<2} {x:>14.8f} {y:>14.8f} {z:>14.8f}"
        for sym, (x, y, z) in zip(atoms.get_chemical_symbols(), atoms.get_positions(), strict=True)
    ]
    # Gaussian requires a trailing blank line after the molecule spec.
    return "\n".join(header + body) + "\n\n"


# ── parsers ────────────────────────────────────────────────────────────────


def is_log_complete(log_path: str) -> bool:
    """True iff log_path exists and ends with a Normal termination marker."""
    if not os.path.exists(log_path):
        return False
    try:
        with open(log_path, "rb") as f:
            f.seek(0, os.SEEK_END)
            size = f.tell()
            f.seek(max(0, size - 4096))
            tail = f.read().decode("utf-8", errors="ignore")
    except OSError:
        return False
    return "Normal termination of Gaussian" in tail


def parse_properties(log_path: str) -> PropertiesResult:
    """Parse HOMO/LUMO (eV), gap, dipole (Debye), SCF energy (eV)."""
    occ, virt = _parse_orbitals_ev(log_path)
    if not occ or not virt:
        raise ValueError(f"{log_path}: no orbital eigenvalues found (need pop=full)")
    homo = occ[-1]
    lumo = virt[0]
    return PropertiesResult(
        homo_ev=homo,
        lumo_ev=lumo,
        gap_ev=lumo - homo,
        dipole_debye=_parse_dipole_debye(log_path),
        energy_ev=_parse_scf_energy_ev(log_path),
    )


def parse_excited_states(log_path: str) -> list[ExcitedState]:
    """Parse TDDFT 'Excited State N' lines into ExcitedState records."""
    states: list[ExcitedState] = []
    with open(log_path) as f:
        for line in f:
            if not line.strip().startswith("Excited State"):
                continue
            parts = line.split()
            try:
                n = int(parts[2].rstrip(":"))
                energy_ev = float(parts[4])
                wavelength_nm = float(parts[6])
                osc = float(next(p for p in parts if p.startswith("f=")).split("=")[1])
            except (IndexError, ValueError, StopIteration):
                continue
            states.append(
                ExcitedState(
                    n=n,
                    energy_ev=energy_ev,
                    wavelength_nm=wavelength_nm,
                    oscillator_strength=osc,
                )
            )
    return sorted(states, key=lambda s: s.energy_ev)


# ── internal log helpers ───────────────────────────────────────────────────


def parse_dimer_orbitals(log_path: str) -> tuple[float, float, float, float]:
    """Return (homo_ev, homo_minus1_ev, lumo_ev, lumo_plus1_ev) from a dimer log.

    Used for the energy-splitting transfer integral:
        J_hole     = |homo_ev - homo_minus1_ev| / 2
        J_electron = |lumo_ev - lumo_plus1_ev|  / 2
    """
    occ, virt = _parse_orbitals_ev(log_path)
    if len(occ) < 2 or len(virt) < 2:
        raise ValueError(f"{log_path}: need ≥2 occ and ≥2 virt eigenvalues (dimer log expected)")
    return occ[-1], occ[-2], virt[0], virt[1]


def _parse_orbitals_ev(log_path: str) -> tuple[list[float], list[float]]:
    occ: list[float] = []
    virt: list[float] = []
    with open(log_path) as f:
        for line in f:
            if "Alpha  occ. eigenvalues" in line:
                occ.extend(float(x) * HARTREE_TO_EV for x in line.split()[4:])
            elif "Alpha virt. eigenvalues" in line:
                virt.extend(float(x) * HARTREE_TO_EV for x in line.split()[4:])
    return occ, virt


def _parse_dipole_debye(log_path: str) -> float | None:
    dipole = None
    in_block = False
    with open(log_path) as f:
        for line in f:
            if "Dipole moment" in line and "Debye" in line:
                in_block = True
            elif in_block:
                if "Tot=" in line:
                    parts = line.split()
                    for i, tok in enumerate(parts):
                        if tok == "Tot=":
                            try:
                                dipole = float(parts[i + 1])
                            except (ValueError, IndexError):
                                pass
                in_block = False
    return dipole


def _parse_scf_energy_ev(log_path: str) -> float:
    last = None
    with open(log_path) as f:
        for line in f:
            if "SCF Done" in line:
                parts = line.split()
                for i, tok in enumerate(parts):
                    if tok == "=":
                        try:
                            last = float(parts[i + 1])
                        except (ValueError, IndexError):
                            pass
                        break
    if last is None:
        raise ValueError(f"{log_path}: no 'SCF Done' line")
    return last * HARTREE_TO_EV
