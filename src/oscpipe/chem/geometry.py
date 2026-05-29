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


def read_gaussian_log(log_path: str):
    """Read the final geometry from a Gaussian .log into an ASE Atoms object."""
    import ase.io

    return ase.io.read(log_path, format="gaussian-out")


def build_pi_stack_dimer(atoms, stack_distance: float = 3.5):
    """Return ASE Atoms: two slip-stacked copies separated by stack_distance Å.

    1. PCA-aligns the monomer to its best-fit π-plane (XY) using heavy atoms only,
       so stack_distance is the true π-plane separation regardless of input orientation.
    2. If cofacial stacking clashes (any inter-monomer distance < 2.5 Å), searches
       a 2D (X, Y) grid sorted by |offset| for the smallest lateral displacement that
       resolves all clashes.  Raises RuntimeError if the grid is exhausted.
    """
    import numpy as np

    # Heavy-atom mask for PCA (exclude H so out-of-plane H doesn't bias the plane normal).
    heavy_mask = np.array([s != "H" for s in atoms.get_chemical_symbols()])
    if heavy_mask.sum() < 3:
        heavy_mask = np.ones(len(atoms), dtype=bool)  # fallback: use all atoms

    pos_all = atoms.positions - atoms.positions.mean(axis=0)
    heavy_pos = pos_all[heavy_mask]
    _, _, Vt = np.linalg.svd(heavy_pos - heavy_pos.mean(axis=0))
    normal = Vt[-1]  # minimum-variance direction = plane normal

    z_hat = np.array([0.0, 0.0, 1.0])
    axis = np.cross(normal, z_hat)
    sin_a = float(np.linalg.norm(axis))
    cos_a = float(np.dot(normal, z_hat))
    if sin_a < 1e-10:
        R = np.eye(3)
    else:
        axis /= sin_a
        K = np.array([[0.0, -axis[2], axis[1]], [axis[2], 0.0, -axis[0]], [-axis[1], axis[0], 0.0]])
        R = np.eye(3) + sin_a * K + (1.0 - cos_a) * (K @ K)

    first = atoms.copy()
    first.positions = (R @ pos_all.T).T
    second = first.copy()
    second.positions[:, 2] += stack_distance

    _MIN_CONTACT = 2.5
    _SLIP_STEP = 0.5
    _MAX_SLIP = 6.5

    # Build grid of (dx, dy) candidates sorted by distance from origin.
    steps = np.arange(-_MAX_SLIP, _MAX_SLIP + _SLIP_STEP * 0.5, _SLIP_STEP)
    candidates = sorted(
        ((dx, dy) for dx in steps for dy in steps),
        key=lambda p: p[0] ** 2 + p[1] ** 2,
    )

    base_pos2 = second.positions.copy()
    for dx, dy in candidates:
        second.positions = base_pos2.copy()
        second.positions[:, 0] += dx
        second.positions[:, 1] += dy
        diffs = first.positions[:, np.newaxis] - second.positions[np.newaxis]
        if np.linalg.norm(diffs, axis=2).min() >= _MIN_CONTACT:
            return first + second

    raise RuntimeError(
        f"build_pi_stack_dimer: could not resolve clashes within ±{_MAX_SLIP} Å grid. "
        "Try a larger stack_distance or manual geometry."
    )
