# ADR 0002 — The dimer SP uses ωB97X-D, not B3LYP (mixed functionals in λ_h)

**Status:** Accepted
**Date:** 2026-05-29
**Authors:** Aarks0rn

## Context

The λ_h workflow (`workflows/lambda_h.py`) runs five Gaussian jobs. Four of them —
the neutral/cation optimisations and the two cross single points — use B3LYP/6-31G\*\*.
The fifth, the π-stacked dimer SP, exists only to extract the hole transfer integral
J_hole, which `analysis.indo.transfer_integral` reads as the HOMO/HOMO-1 eigenvalue
splitting (`|E_HOMO − E_{HOMO-1}| / 2`) at a **fixed** geometry.

The workspace `CLAUDE.md` previously listed two open dimer issues:

1. *Clashing monomer geometries → unphysical "+40 eV binding energies."*
2. *Dimer SP has no dispersion — switch to B3LYP-D3 or ωB97X-D.*

On inspection both descriptions predate the current code:

- The clash bug is closed. `chem.geometry.build_pi_stack_dimer` PCA-aligns the monomer
  to its π-plane and slip-stacks on an (x, y) grid until no inter-monomer contact is
  below 2.5 Å (tests `test_build_pi_stack_dimer_no_clash_*`). No binding/interaction
  energy is computed anywhere in the pipeline — the dimer log is consumed solely for
  its eigenvalues — so the "+40 eV" figure refers to a computation that no longer exists.
- "Add dispersion" as literally stated would be a **no-op for J**. Grimme D3
  (`EmpiricalDispersion=GD3`) is an additive, atom-pairwise correction applied after SCF
  convergence; it depends only on nuclear positions, never enters the Kohn-Sham matrix,
  and therefore leaves the orbital eigenvalues — and thus the splitting J is read from —
  unchanged. The only lever on J accuracy at a fixed geometry is the **functional itself**.

## Decision

Run the dimer SP with **ωB97X-D/6-31G\*\*** while the rest of the λ_h workflow stays on
**B3LYP/6-31G\*\***.

ωB97X-D is a range-separated hybrid: its exact long-range exchange reduces the
self-interaction / delocalization error that inflates B3LYP HOMO/HOMO-1 splittings for
π-stacked dimers. The monomer copies are identical (symmetric homodimer), so the
energy-splitting identity J = ΔE/2 is valid and the functional is the only thing that
changes J.

We deliberately do **not** switch the four λ_h jobs to ωB97X-D. Doing so would change
the reorganisation-energy reference values (and break the 2.72114 eV regression test) for
no benefit to this fix; the reorganisation-energy scheme is a separate accuracy question.

## Consequences

- One workflow now spans two functionals. This is intentional and documented here so a
  future reviewer does not "fix the inconsistency" by collapsing them. The job signature
  (`store.cache.signature`) includes the method, so the dimer job caches independently of
  the B3LYP jobs — no cache collision.
- λ_h, the four-point energies, and the test suite are unchanged (112 pass / 2 skip); the
  test double returns synthetic logs by label suffix, independent of the method string.
- **Do not re-add a D3 dispersion correction to the dimer SP expecting it to change J** —
  it cannot, for the reason in Context. If a future workflow needs an actual π-stacking
  *binding* energy (not J), that is a new computation and warrants its own ADR.

## Alternatives rejected

- **Keep B3LYP, add `EmpiricalDispersion=GD3`.** No-op for J (see Context). Rejected.
- **Switch the whole workflow to ωB97X-D.** Consistent, but re-baselines λ_h and breaks
  the regression test for a property the fix was not about. Rejected — out of scope.

## References

- `workflows/lambda_h.py` — dimer SP block (`dimer_method = "wb97xd"`)
- `analysis/indo.py`, `dft.gaussian.parse_dimer_orbitals` — J from eigenvalue splitting
- `chem/geometry.py` — `build_pi_stack_dimer` (clash-free stacking)
- ADR 0001 — orchestration layers (where this workflow lives)
