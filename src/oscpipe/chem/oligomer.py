"""Stitch a [*]…[*] repeat-unit SMILES into a linear n-mer.

A conjugated-polymer repeat unit is written as a SMILES with exactly two ``[*]``
dummy atoms marking the backbone attachment points; in SMILES source order they
are (head, tail). ``build_oligomer`` concatenates ``n`` copies tail->head and
caps the terminal head/tail with H, returning the canonical n-mer SMILES.

This is the geometry input for the oligomer-length sweep: run n = 1, 2, 3 of an
anchor's repeat unit, then extrapolate a property vs 1/n to the polymer limit
(see ``oscpipe.analysis.extrapolation``). Pure function — no IO, no DFT.
"""

from __future__ import annotations

from rdkit import Chem


def build_oligomer(repeat_unit_smiles: str, n: int) -> str:
    """Return the canonical SMILES of the n-mer of ``repeat_unit_smiles``.

    ``repeat_unit_smiles`` must contain exactly two ``[*]`` attachment points.
    Raises ``ValueError`` on an unparseable SMILES, the wrong number of
    attachment points, or ``n < 1``.
    """
    if n < 1:
        raise ValueError(f"n must be >= 1, got {n}")
    unit = Chem.MolFromSmiles(repeat_unit_smiles)
    if unit is None:
        raise ValueError(f"could not parse repeat-unit SMILES: {repeat_unit_smiles!r}")

    dummies = [a.GetIdx() for a in unit.GetAtoms() if a.GetAtomicNum() == 0]
    if len(dummies) != 2:
        raise ValueError(
            f"repeat unit needs exactly two [*] attachment points, found {len(dummies)}"
        )
    head_dummy, tail_dummy = dummies
    head_nbr = unit.GetAtomWithIdx(head_dummy).GetNeighbors()[0].GetIdx()
    tail_nbr = unit.GetAtomWithIdx(tail_dummy).GetNeighbors()[0].GetIdx()

    width = unit.GetNumAtoms()
    chain = unit
    for _ in range(n - 1):
        chain = Chem.CombineMols(chain, unit)

    rw = Chem.RWMol(chain)
    # Join unit i's tail-neighbour to unit (i+1)'s head-neighbour.
    for i in range(n - 1):
        rw.AddBond(i * width + tail_nbr, (i + 1) * width + head_nbr, Chem.BondType.SINGLE)
    # Drop every dummy: internal ones are replaced by the new bond, terminal
    # ones become implicit-H caps. Remove high indices first to keep them valid.
    drop = [i * width + d for i in range(n) for d in (head_dummy, tail_dummy)]
    for idx in sorted(drop, reverse=True):
        rw.RemoveAtom(idx)

    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    return Chem.MolToSmiles(mol)
