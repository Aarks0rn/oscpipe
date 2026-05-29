"""λ_h workflow — hole reorganisation energy + transfer integral + Marcus rate.

Composes the five Gaussian jobs (4-point Nelsen scheme + π-stacked dimer SP) via
the :class:`oscpipe.run.JobRunner`, then derives the physics with the tested
``analysis`` functions.  The arithmetic lives in ``analysis.marcus`` /
``analysis.indo``, not inline here — this module only sequences the jobs.

Layer: above ``run`` / ``analysis``; the CLI's ``lambda`` subcommand is a thin
adapter onto :func:`run_lambda_h`.
"""

from __future__ import annotations

import json
import sys

from ..analysis.indo import transfer_integral
from ..analysis.marcus import lambda_hole_from_4_points, marcus_rate
from ..chem import smiles as chem_smiles
from ..chem.geometry import build_pi_stack_dimer, read_gaussian_log
from ..dft import gaussian
from ..run import JobRunner, _now
from ..run import label as _label
from ..settings import Settings
from ..store import db
from ..store.cache import signature


def run_lambda_h(
    args,
    settings: Settings,
    backend,
    conn,
    *,
    stdout=None,
    geometry_loader=None,
) -> int:
    """Submit the 4 Gaussian jobs for λ_h (hole reorganisation energy).

    Sequence:
        1. neutral_opt  (charge=0, mult=1) → E_N(N) and neutral-opt geom
        2. cation_opt   (charge=+1, mult=2) → E_C(C) and cation-opt geom
        3. neutral_sp at cation geom → E_N(C)
        4. cation_sp at neutral geom → E_C(N)

    λ_h = [E_C(N) − E_C(C)] + [E_N(C) − E_N(N)]
    """
    if stdout is None:
        stdout = sys.stdout
    if geometry_loader is None:
        geometry_loader = read_gaussian_log

    canonical, warnings = chem_smiles.canonicalise(args.smiles)
    for w in warnings:
        print(f"warning: {w}", file=stdout)

    initial = chem_smiles.embed_3d(canonical)
    method, basis = "b3lyp", "6-31g**"

    wf_id = db.insert_workflow(conn, "lambda_h", canonical, _now())
    print(f"workflow: id={wf_id} kind=lambda_h smiles={canonical}", file=stdout)

    runner = JobRunner(conn, settings, backend, stdout)

    def _opt_job(charge, mult, suffix):
        sig = signature(canonical, method, basis, charge, mult)
        label = f"{_label(canonical, sig)}_{suffix}"

        def _build():
            return gaussian.write_com_properties(
                initial,
                method=method,
                basis=basis,
                charge=charge,
                mult=mult,
                nproc=settings.gaussian_nproc,
                mem=settings.gaussian_mem,
                label=label,
                chk=f"{label}.chk",
            )

        job = db.Job(
            id=None,
            signature=sig,
            smiles=canonical,
            method=method,
            basis=basis,
            charge=charge,
            mult=mult,
            job_kind="properties",
            status="pending",
            submitted_at=_now(),
            workflow_id=wf_id,
            notes=suffix,
        )
        # need_log: the optimised geometry is read back from the log for the SPs.
        return runner.resolve(job, label, _build, wait=True, need_log=True)

    # 1+2: optimise neutral and cation.
    neutral = _opt_job(0, 1, "neutral_opt")
    cation = _opt_job(1, 2, "cation_opt")
    if neutral.status != "complete" or cation.status != "complete":
        db.update_workflow(
            conn,
            wf_id,
            "error",
            summary_json=json.dumps(
                {"stage": "opt", "neutral": neutral.status, "cation": cation.status}
            ),
        )
        print("lambda: opt stage failed", file=stdout)
        return 1
    e_neut_opt = neutral.result.energy_ev
    e_cat_opt = cation.result.energy_ev
    neutral_geom = geometry_loader(neutral.log_path)
    cation_geom = geometry_loader(cation.log_path)

    # 3+4: single points at the *other* geometry.
    def _sp_job(atoms, charge, mult, kind, suffix):
        sig = signature(canonical, method, basis, charge, mult, extras=suffix)
        label = f"{_label(canonical, sig)}_{suffix}"

        def _build():
            return gaussian.write_com_sp(
                atoms,
                method=method,
                basis=basis,
                charge=charge,
                mult=mult,
                nproc=settings.gaussian_nproc,
                mem=settings.gaussian_mem,
                label=label,
                chk=f"{label}.chk",
            )

        job = db.Job(
            id=None,
            signature=sig,
            smiles=canonical,
            method=method,
            basis=basis,
            charge=charge,
            mult=mult,
            job_kind=kind,
            status="pending",
            submitted_at=_now(),
            workflow_id=wf_id,
            notes=suffix,
        )
        # The SP contributes only its energy (read from the result), not geometry.
        return runner.resolve(job, label, _build, wait=True)

    n_at_c = _sp_job(cation_geom, 0, 1, "sp_neutral", "neutral_at_cation_geom")
    c_at_n = _sp_job(neutral_geom, 1, 2, "sp_cation", "cation_at_neutral_geom")
    if n_at_c.status != "complete" or c_at_n.status != "complete":
        db.update_workflow(
            conn,
            wf_id,
            "error",
            summary_json=json.dumps(
                {"stage": "sp", "n_at_c": n_at_c.status, "c_at_n": c_at_n.status}
            ),
        )
        print("lambda: sp stage failed", file=stdout)
        return 1

    e_neut_at_cat = n_at_c.result.energy_ev
    e_cat_at_neut = c_at_n.result.energy_ev

    # 4-point Nelsen scheme (energies already in eV). E_n_c is the cation-charge
    # SP at the neutral geometry; E_c_n is the neutral-charge SP at the cation
    # geometry — see analysis.marcus for the mapping.
    lambda_h = lambda_hole_from_4_points(
        e_n_n=e_neut_opt,
        e_n_c=e_cat_at_neut,
        e_c_c=e_cat_opt,
        e_c_n=e_neut_at_cat,
    )
    print(f"lambda_h = {lambda_h:.4f} eV", file=stdout)

    # ── job 5: dimer SP → J_hole → Marcus rate ──────────────────────────────
    # The dimer SP uses ωB97X-D, not the B3LYP used for the λ_h jobs above. J_hole
    # is the HOMO/HOMO-1 eigenvalue splitting at a fixed geometry, so the functional
    # is the only lever on its accuracy: ωB97X-D's range separation reduces the
    # delocalization error that inflates B3LYP splittings. (A B3LYP-D3 dispersion
    # correction would be a no-op here — D3 is a post-SCF additive term that never
    # enters the Kohn-Sham matrix, so it leaves the eigenvalues unchanged.) See
    # docs/adr/0002-dimer-functional.md for the mixed-functional rationale.
    dimer_method = "wb97xd"
    dimer_atoms = build_pi_stack_dimer(neutral_geom)
    sig_dimer = signature(canonical, dimer_method, basis, 0, 1, extras="dimer_sp")
    label_dimer = f"{_label(canonical, sig_dimer)}_dimer_sp"

    def _build_dimer():
        return gaussian.write_com_sp(
            dimer_atoms,
            method=dimer_method,
            basis=basis,
            charge=0,
            mult=1,
            nproc=settings.gaussian_nproc,
            mem=settings.gaussian_mem,
            label=label_dimer,
            chk=f"{label_dimer}.chk",
        )

    job_dimer = db.Job(
        id=None,
        signature=sig_dimer,
        smiles=canonical,
        method=dimer_method,
        basis=basis,
        charge=0,
        mult=1,
        job_kind="sp_dimer",
        status="pending",
        submitted_at=_now(),
        workflow_id=wf_id,
        notes="dimer_sp",
    )
    # need_log: J_hole is read from the dimer log by analysis.indo.transfer_integral.
    dimer = runner.resolve(job_dimer, label_dimer, _build_dimer, wait=True, need_log=True)

    j_hole_ev: float | None = None
    marcus: float | None = None
    if dimer.status == "complete" and dimer.log_path:
        j_hole_ev = transfer_integral(dimer.log_path)
        marcus = marcus_rate(lambda_h, j_hole_ev, delta_g_ev=0.0)
        print(f"J_hole = {j_hole_ev:.4f} eV", file=stdout)
        print(f"marcus_rate = {marcus:.4e} s^-1", file=stdout)
    else:
        print(f"dimer SP {dimer.status} — J/marcus not available", file=stdout)

    summary = {
        "lambda_h_ev": lambda_h,
        "j_hole_ev": j_hole_ev,
        "marcus_rate_s1": marcus,
        "energies_ev": {
            "neutral_opt": e_neut_opt,
            "cation_opt": e_cat_opt,
            "neutral_at_cation_geom": e_neut_at_cat,
            "cation_at_neutral_geom": e_cat_at_neut,
        },
    }
    db.update_workflow(conn, wf_id, "complete", summary_json=json.dumps(summary))
    print(f"workflow id={wf_id}", file=stdout)
    return 0
