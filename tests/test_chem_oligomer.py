"""Tests for chem/oligomer: stitch a [*]…[*] repeat unit into an n-mer.

Pure-function, fixture-free. RDKit runs locally; no DFT. Convention under test:
the two [*] dummy atoms, in SMILES source order, are (head, tail); units are
joined tail->head and the terminal head/tail are capped with H.
"""

import pytest
from rdkit import Chem

from oscpipe.chem.oligomer import build_oligomer

THIOPHENE_25 = "[*]c1ccc([*])s1"  # 2,5-thiophene-diyl repeat unit


def _canon(s: str) -> str:
    return Chem.CanonSmiles(s)


def test_n1_caps_both_attachment_points_to_bare_monomer():
    # n=1: both [*] become H, leaving plain thiophene.
    assert _canon(build_oligomer(THIOPHENE_25, 1)) == _canon("c1ccsc1")


def test_n2_links_two_units_into_bithiophene():
    assert _canon(build_oligomer(THIOPHENE_25, 2)) == _canon("c1ccc(-c2cccs2)s1")


def test_n3_has_three_units_worth_of_heavy_atoms_and_no_dummies():
    mol = Chem.MolFromSmiles(build_oligomer(THIOPHENE_25, 3))
    assert mol is not None
    assert sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 0) == 0
    assert mol.GetNumAtoms() == 3 * 5  # C4S per thiophene ring


def test_asymmetric_unit_n2_preserves_formula():
    # 3-fluoro-2,5-thiophene-diyl: two units, two F, no dummies.
    mol = Chem.MolFromSmiles(build_oligomer("[*]c1cc(F)c([*])s1", 2))
    assert mol is not None
    counts = {}
    for a in mol.GetAtoms():
        counts[a.GetSymbol()] = counts.get(a.GetSymbol(), 0) + 1
    assert counts.get("*", 0) == 0
    assert counts["F"] == 2
    assert counts["S"] == 2
    assert counts["C"] == 8


def test_requires_exactly_two_attachment_points():
    with pytest.raises(ValueError):
        build_oligomer("c1ccccc1", 2)  # zero [*]


def test_rejects_nonpositive_n():
    with pytest.raises(ValueError):
        build_oligomer(THIOPHENE_25, 0)


def test_invalid_smiles_raises():
    with pytest.raises(ValueError):
        build_oligomer("definitely not smiles", 1)
