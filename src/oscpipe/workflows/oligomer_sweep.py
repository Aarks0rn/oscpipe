"""Oligomer-length sweep — extrapolate HOMO / gap / optical gap to the polymer limit.

For each repeat-unit count n = 1..max_n, build the n-mer SMILES
(:func:`oscpipe.chem.oligomer.build_oligomer`), run a B3LYP ``properties`` opt
(HOMO + KS gap) and an ωB97X-D ``tddft`` job at that optimised geometry (optical
gap = lowest excited state), then extrapolate each property vs 1/n to the
n -> infinity limit (:func:`oscpipe.analysis.extrapolation.extrapolate_inverse_n`).

Cheap tier only — no λ_h / dimer here; those run later on screen survivors. The
``properties`` job is ``opt pop=full`` with no frequency, which keeps the
large-oligomer geometry optimisations affordable. Layer: above ``run`` /
``analysis`` / ``chem``; the CLI ``oligomer`` subcommand is a thin adapter onto
:func:`run_oligomer_sweep`.
"""

from __future__ import annotations

import json
import sys

from ..analysis.extrapolation import extrapolate_inverse_n
from ..chem import smiles as chem_smiles
from ..chem.geometry import read_gaussian_log
from ..chem.oligomer import build_oligomer
from ..dft import gaussian
from ..run import JobRunner, _now
from ..run import label as _label
from ..settings import Settings
from ..store import db
from ..store.cache import signature

PREOPT_METHOD = "pm6"  # cheap semi-empirical geometry seed (no basis)
OPT_METHOD = "b3lyp"
OPTICAL_METHOD = "wb97xd"
BASIS = "6-31g**"
_EXTRAPOLATED_KEYS = ("homo_ev", "gap_ks_ev", "optical_gap_ev")


def run_oligomer_sweep(
    args,
    settings: Settings,
    backend,
    conn,
    *,
    stdout=None,
    geometry_loader=None,
) -> int:
    """Run the n = 1..max_n sweep and extrapolate each property to the polymer limit."""
    if stdout is None:
        stdout = sys.stdout
    if geometry_loader is None:
        geometry_loader = read_gaussian_log

    repeat_unit, max_n, nstates = args.repeat_unit, args.max_n, args.nstates
    if max_n < 2:
        print(f"oligomer: need max_n >= 2 to extrapolate, got {max_n}", file=stdout)
        return 1

    wf_id = db.insert_workflow(conn, "oligomer_sweep", repeat_unit, _now())
    print(
        f"workflow: id={wf_id} kind=oligomer_sweep repeat_unit={repeat_unit} max_n={max_n}",
        file=stdout,
    )
    runner = JobRunner(conn, settings, backend, stdout)

    per_n: list[dict] = []
    for n in range(1, max_n + 1):
        canonical, warnings = chem_smiles.canonicalise(build_oligomer(repeat_unit, n))
        for w in warnings:
            print(f"warning (n={n}): {w}", file=stdout)
        initial = chem_smiles.embed_3d(canonical)

        # 0. PM6 pre-opt — seed a near-physical geometry. A direct B3LYP/6-31G**
        #    opt on a large n=3 oligomer thrashes (huge steps, ~14 h/step, never
        #    converging); pre-optimising at PM6 first cuts the B3LYP opt to a few
        #    steps. The optimised PM6 geometry replaces the RDKit embed below.
        sig_pre = signature(canonical, PREOPT_METHOD, "", 0, 1, job_kind="preopt")
        label_pre = f"{_label(canonical, sig_pre)}_n{n}_preopt"

        def _build_preopt(atoms=initial, lbl=label_pre):
            return gaussian.write_com_preopt(
                atoms,
                charge=0,
                mult=1,
                nproc=settings.gaussian_nproc,
                mem=settings.gaussian_mem,
                label=lbl,
                chk=f"{lbl}.chk",
            )

        job_pre = db.Job(
            id=None,
            signature=sig_pre,
            smiles=canonical,
            method=PREOPT_METHOD,
            basis="",
            charge=0,
            mult=1,
            job_kind="preopt",
            status="pending",
            submitted_at=_now(),
            workflow_id=wf_id,
            notes=f"n{n}_preopt",
        )
        pre = runner.resolve(job_pre, label_pre, _build_preopt, wait=True, need_log=True)
        if pre.status != "complete":
            return _fail(conn, wf_id, stdout, "preopt", n, pre.status)
        initial = geometry_loader(pre.log_path)  # B3LYP opt starts from the PM6 geometry

        # 1. properties (B3LYP opt) — HOMO, KS gap, optimised geometry.
        sig_p = signature(canonical, OPT_METHOD, BASIS, 0, 1)
        label_p = f"{_label(canonical, sig_p)}_n{n}_opt"

        def _build_props(atoms=initial, lbl=label_p):
            return gaussian.write_com_properties(
                atoms,
                method=OPT_METHOD,
                basis=BASIS,
                charge=0,
                mult=1,
                nproc=settings.gaussian_nproc,
                mem=settings.gaussian_mem,
                label=lbl,
                chk=f"{lbl}.chk",
                # Floppy TVT backbones thrash a default-trust-radius opt (job 71248);
                # damp it and accept loose convergence — HOMO/gap/optical are
                # insensitive to the last bit of geometry. See write_com_properties.
                opt_route="opt=(loose,maxstep=10,maxcycles=300)",
            )

        job_p = db.Job(
            id=None,
            signature=sig_p,
            smiles=canonical,
            method=OPT_METHOD,
            basis=BASIS,
            charge=0,
            mult=1,
            job_kind="properties",
            status="pending",
            submitted_at=_now(),
            workflow_id=wf_id,
            notes=f"n{n}_opt",
        )
        prop = runner.resolve(job_p, label_p, _build_props, wait=True, need_log=True)
        if prop.status != "complete":
            return _fail(conn, wf_id, stdout, "opt", n, prop.status)
        opt_geom = geometry_loader(prop.log_path)

        # 2. tddft (ωB97X-D) at the optimised geometry — optical gap = lowest state.
        sig_t = signature(canonical, OPTICAL_METHOD, BASIS, 0, 1, extras=f"tddft{nstates}")
        label_t = f"{_label(canonical, sig_t)}_n{n}_tddft"

        def _build_tddft(atoms=opt_geom, lbl=label_t):
            return gaussian.write_com_tddft(
                atoms,
                method=OPTICAL_METHOD,
                basis=BASIS,
                charge=0,
                mult=1,
                nstates=nstates,
                nproc=settings.gaussian_nproc,
                mem=settings.gaussian_mem,
                label=lbl,
                chk=f"{lbl}.chk",
            )

        job_t = db.Job(
            id=None,
            signature=sig_t,
            smiles=canonical,
            method=OPTICAL_METHOD,
            basis=BASIS,
            charge=0,
            mult=1,
            job_kind="tddft",
            status="pending",
            submitted_at=_now(),
            workflow_id=wf_id,
            notes=f"n{n}_tddft",
        )
        tdd = runner.resolve(job_t, label_t, _build_tddft, wait=True, need_log=True)
        if tdd.status != "complete" or not tdd.log_path:
            return _fail(conn, wf_id, stdout, "tddft", n, tdd.status)
        states = gaussian.parse_excited_states(tdd.log_path)
        if not states:
            return _fail(conn, wf_id, stdout, "tddft", n, "no excited states")
        optical_gap = states[0].energy_ev  # parse_excited_states sorts by energy → S1

        per_n.append(
            {
                "n": n,
                "smiles": canonical,
                "homo_ev": prop.result.homo_ev,
                "gap_ks_ev": prop.result.gap_ev,
                "optical_gap_ev": optical_gap,
            }
        )
        print(
            f"  n={n}: HOMO={prop.result.homo_ev:.3f} gap_ks={prop.result.gap_ev:.3f} "
            f"optical={optical_gap:.3f} eV",
            file=stdout,
        )

    extrapolated = {}
    for key in _EXTRAPOLATED_KEYS:
        r = extrapolate_inverse_n([(d["n"], d[key]) for d in per_n])
        extrapolated[key] = {"limit": r.limit, "slope": r.slope, "r_squared": r.r_squared}
        print(f"  {key} limit = {r.limit:.3f} eV (R^2={r.r_squared:.4f})", file=stdout)

    summary = {
        "repeat_unit": repeat_unit,
        "max_n": max_n,
        "nstates": nstates,
        "per_n": per_n,
        "extrapolated": extrapolated,
    }
    db.update_workflow(conn, wf_id, "complete", summary_json=json.dumps(summary))
    print(f"workflow id={wf_id}", file=stdout)
    return 0


def _fail(conn, wf_id: int, stdout, stage: str, n: int, status: str) -> int:
    db.update_workflow(
        conn, wf_id, "error", summary_json=json.dumps({"stage": stage, "n": n, "status": status})
    )
    print(f"oligomer: {stage} failed at n={n} ({status})", file=stdout)
    return 1
