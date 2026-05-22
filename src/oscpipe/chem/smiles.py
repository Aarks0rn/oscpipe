"""SMILES sanitisation + 3D embedding.

Public API:
    canonicalise(smiles) -> (canonical_smiles, warnings)
    embed_3d(smiles) -> ase.Atoms

RDKit is the primary path; openbabel is the fallback for SMILES that RDKit
can parse-but-not-embed. The RDKit path matters for fused polycyclic
aromatics where openbabel's make3D puckers the ring (the ChDT bug).
"""

from __future__ import annotations

import tempfile
from pathlib import Path


def canonicalise(smiles: str) -> tuple[str, list[str]]:
    """Canonicalise a SMILES via RDKit.

    Returns ``(canonical_smiles, warnings)``. On parse failure the original
    string is returned with a warning so callers may still try openbabel.
    Warnings are also emitted for sp3 carbons sandwiched between aromatic
    neighbours (the ChDT case — almost always a SMILES-authoring mistake).
    """
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return smiles, [f"RDKit could not parse SMILES {smiles!r}; geometry may be wrong."]

    warnings: list[str] = []

    sp3_in_arom_ring = [
        atom.GetIdx()
        for atom in mol.GetAtoms()
        if atom.GetSymbol() == "C"
        and atom.GetHybridization() == Chem.HybridizationType.SP3
        and atom.IsInRing()
        and sum(1 for nb in atom.GetNeighbors() if nb.GetIsAromatic()) >= 2
    ]
    if sp3_in_arom_ring:
        warnings.append(
            f"sp3 carbon(s) at atom index {sp3_in_arom_ring} sit between "
            f"aromatic neighbours; DFT optimisation will pucker them out of "
            f"plane. If the molecule is meant to be planar, redraw with "
            f"aromatic bonds (lowercase SMILES)."
        )

    canonical = Chem.MolToSmiles(mol)
    if canonical != smiles:
        warnings.append(f"SMILES canonicalised: {smiles!r} → {canonical!r}.")
    return canonical, warnings


def embed_3d(smiles: str):
    """Return an ``ase.Atoms`` with 3D coordinates for ``smiles``.

    Tries RDKit (``EmbedMolecule`` + ``MMFFOptimizeMolecule``) first, falling
    back to openbabel's ``make3D`` if RDKit cannot embed. Raises ``ValueError``
    if neither path produces a geometry.
    """
    import ase.io

    with tempfile.TemporaryDirectory() as tmp:
        sdf = Path(tmp) / "embed.sdf"
        if _rdkit_embed_to_sdf(smiles, str(sdf)):
            return ase.io.read(str(sdf))
        if _openbabel_embed_to_sdf(smiles, str(sdf)):
            return ase.io.read(str(sdf))
    raise ValueError(f"could not embed SMILES {smiles!r} via RDKit or openbabel")


def _rdkit_embed_to_sdf(smiles: str, sdf_file: str, max_iters: int = 500) -> bool:
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
        return False
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=max_iters)
    except Exception:
        try:
            AllChem.UFFOptimizeMolecule(mol, maxIters=max_iters)
        except Exception:
            return False
    with Chem.SDWriter(sdf_file) as w:
        w.write(mol)
    return True


def _openbabel_embed_to_sdf(smiles: str, sdf_file: str) -> bool:
    try:
        from openbabel import pybel
    except ImportError:
        return False
    try:
        mol = pybel.readstring("smi", smiles)
    except OSError:
        return False
    mol.make3D()
    mol.write("sdf", sdf_file, overwrite=True)
    return True
