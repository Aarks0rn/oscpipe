"""Tests for `oscpipe.chem.geometry` atoms ↔ xyz helpers."""

from __future__ import annotations

import ase
import numpy as np

from oscpipe.chem.geometry import atoms_to_xyz, xyz_to_atoms


def test_atoms_to_xyz_roundtrip_preserves_geometry():
    atoms = ase.Atoms("H2O", positions=[(0, 0, 0), (0.96, 0, 0), (-0.24, 0.93, 0)])
    xyz = atoms_to_xyz(atoms)
    restored = xyz_to_atoms(xyz)
    assert restored.get_chemical_symbols() == [
        "H",
        "H",
        "O",
    ] or restored.get_chemical_symbols() == list(atoms.get_chemical_symbols())
    np.testing.assert_allclose(restored.get_positions(), atoms.get_positions(), atol=1e-6)


def test_atoms_to_xyz_header_lines():
    atoms = ase.Atoms("H2", positions=[(0, 0, 0), (0, 0, 0.74)])
    xyz = atoms_to_xyz(atoms)
    lines = xyz.splitlines()
    assert lines[0] == "2"
    assert lines[1] == ""
    assert lines[2].startswith("H ")
    assert lines[3].startswith("H ")


def test_xyz_to_atoms_parses_with_comment():
    xyz = "2\ncomment line\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n"
    atoms = xyz_to_atoms(xyz)
    assert atoms.get_chemical_symbols() == ["H", "H"]
    np.testing.assert_allclose(atoms.get_positions()[1], [0.0, 0.0, 0.74], atol=1e-6)


def test_geometry_hash_stable_and_distinct():
    from oscpipe.chem.geometry import geometry_hash

    a = ase.Atoms(symbols=["C", "H"], positions=[[0, 0, 0], [1.1, 0, 0]])
    b = ase.Atoms(symbols=["C", "H"], positions=[[0, 0, 0], [1.1, 0, 0]])
    c = ase.Atoms(symbols=["C", "H"], positions=[[0, 0, 0], [1.4, 0, 0]])
    assert geometry_hash(a) == geometry_hash(b)
    assert geometry_hash(a) != geometry_hash(c)
    assert len(geometry_hash(a)) == 8
