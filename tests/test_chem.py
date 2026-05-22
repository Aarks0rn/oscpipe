"""Tests for chem/smiles canonicalisation + 3D embed.

Pure-function fixture-free tests. RDKit/openbabel run locally; no DFT.
"""

import numpy as np
import pytest

from oscpipe.chem import smiles

# ── canonicalise ───────────────────────────────────────────────────────────


def test_canonicalise_benzene_lowercase():
    canon, warnings = smiles.canonicalise("c1ccccc1")
    assert canon == "c1ccccc1"
    assert warnings == []


def test_canonicalise_kekulised_form_to_aromatic():
    # Kekulé form should canonicalise to the aromatic form; this is informational.
    canon, warnings = smiles.canonicalise("C1=CC=CC=C1")
    assert canon == "c1ccccc1"
    # First warning announces the rewrite.
    assert any("canonicalised" in w for w in warnings)


def test_canonicalise_flags_sp3_in_aromatic_ring():
    # Fluorene: central methylene is sp3 in a ring, with two aromatic
    # neighbours. Drawing it as a fully aromatic system is a common
    # SMILES-authoring mistake on planar π systems.
    fluorene = "c1ccc2c(c1)Cc1ccccc1-2"
    _, warnings = smiles.canonicalise(fluorene)
    assert any("sp3 carbon" in w for w in warnings)


def test_canonicalise_invalid_returns_original():
    bad = "this is not a smiles"
    canon, warnings = smiles.canonicalise(bad)
    assert canon == bad
    assert any("could not parse" in w.lower() for w in warnings)


# ── embed_3d ───────────────────────────────────────────────────────────────


def test_embed_3d_benzene_atom_count():
    atoms = smiles.embed_3d("c1ccccc1")
    # Benzene: 6 C + 6 H after AddHs
    assert len(atoms) == 12
    symbols = sorted(atoms.get_chemical_symbols())
    assert symbols.count("C") == 6
    assert symbols.count("H") == 6


def test_embed_3d_benzene_planar():
    atoms = smiles.embed_3d("c1ccccc1")
    # Ring carbons should sit close to a plane. Take the 6 carbons by symbol,
    # fit a plane via SVD, residual should be small.
    c_positions = np.array(
        [
            p
            for p, s in zip(atoms.get_positions(), atoms.get_chemical_symbols(), strict=True)
            if s == "C"
        ]
    )
    centred = c_positions - c_positions.mean(axis=0)
    _, sv, _ = np.linalg.svd(centred, full_matrices=False)
    # Smallest singular value ~ out-of-plane deviation.
    assert sv[-1] < 0.3, f"ring is too puckered: sv={sv}"


def test_embed_3d_invalid_raises():
    with pytest.raises(ValueError):
        smiles.embed_3d("this is not a smiles")
