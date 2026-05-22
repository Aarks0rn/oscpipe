"""ASE Atoms ↔ XYZ helpers."""

from __future__ import annotations


def atoms_to_xyz(atoms) -> str:
    """Serialise ASE Atoms to a standard XYZ-format string.

    Line 1 = atom count, line 2 = comment, then `symbol x y z` per atom.
    """
    symbols = atoms.get_chemical_symbols()
    positions = atoms.get_positions()
    lines = [str(len(symbols)), ""]
    for sym, (x, y, z) in zip(symbols, positions):
        lines.append(f"{sym} {x:.8f} {y:.8f} {z:.8f}")
    return "\n".join(lines) + "\n"


def xyz_to_atoms(xyz: str):
    """Parse XYZ-format string into an ASE Atoms object."""
    import ase

    lines = xyz.strip().splitlines()
    n = int(lines[0].strip())
    symbols: list[str] = []
    positions: list[tuple[float, float, float]] = []
    for line in lines[2 : 2 + n]:
        parts = line.split()
        symbols.append(parts[0])
        positions.append((float(parts[1]), float(parts[2]), float(parts[3])))
    return ase.Atoms(symbols=symbols, positions=positions)
