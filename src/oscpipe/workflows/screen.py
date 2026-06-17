"""Screen workflow — the full pipeline over a list of candidates.

Per candidate, sequentially and blocking: properties opt → TD-DFT → the λ_h
workflow.  This module only composes existing pieces (``JobRunner.resolve`` and
``run_lambda_h``) — no new physics.  A failed candidate is recorded and the
campaign continues with the next one.

Resume = re-run the same list: every completed sub-job is a signature-cache hit
(ADR 0003), so only the missing work is submitted.  There is no state machine.
See docs/adr/0004-screen-campaign-workflow.md.

Layer: above ``run``, beside ``lambda_h``; the CLI's ``screen`` subcommand is a
thin adapter onto :func:`run_screen`.
"""

from __future__ import annotations

import json
import sys
from types import SimpleNamespace

from ..chem import smiles as chem_smiles
from ..dft import gaussian
from ..run import JobRunner, _now, make_job
from ..settings import Settings
from ..store import db
from .lambda_h import run_lambda_h


def _mark(r) -> str:
    """Candidate-progress tag for one Resolved: cached / ok / the raw status."""
    if r.cached:
        return "cached"
    return "ok" if r.status == "complete" else r.status


def run_screen(
    smiles_list: list[str],
    args,
    settings: Settings,
    backend,
    conn,
    *,
    stdout=None,
    geometry_loader=None,
) -> int:
    """Run properties → TD-DFT → λ_h for each SMILES in ``smiles_list``.

    ``args`` carries ``method``, ``basis`` and ``nstates`` (the CLI namespace).
    Returns 0 when every candidate completed, 1 otherwise.
    """
    if stdout is None:
        stdout = sys.stdout

    wf_id = db.insert_workflow(
        conn, "screen", f"screen:{len(smiles_list)}", _now(), status="running"
    )
    print(f"screen: workflow_id={wf_id} {len(smiles_list)} candidates", file=stdout)

    runner = JobRunner(conn, settings, backend, stdout)
    progress: list[dict] = []
    failures = 0

    for i, smi in enumerate(smiles_list, 1):
        # Appended before any work so a crash mid-candidate still leaves a record.
        entry: dict = {"smiles": smi, "properties": None, "tddft": None, "lambda_h": None}
        progress.append(entry)
        print(f"screen: [{i}/{len(smiles_list)}] {smi}", file=stdout)
        try:
            canonical, warnings = chem_smiles.canonicalise(smi)
            for w in warnings:
                print(f"warning: {w}", file=stdout)
            entry["smiles"] = canonical

            # Step 1 — properties opt. Same signature form as `oscpipe submit`;
            # at the default method/basis it also equals λ_h's neutral_opt
            # signature, so step 3 cache-hits this job instead of re-optimising.
            job, label = make_job(
                canonical, args.method, args.basis, 0, 1, "properties", workflow_id=wf_id
            )

            def _build_props(canonical=canonical, label=label):
                atoms = chem_smiles.embed_3d(canonical)
                return gaussian.write_com_properties(
                    atoms,
                    method=args.method,
                    basis=args.basis,
                    charge=0,
                    mult=1,
                    nproc=settings.gaussian_nproc,
                    mem=settings.gaussian_mem,
                    label=label,
                    chk=f"{label}.chk",
                )

            r = runner.resolve(job, label, _build_props, wait=True)
            entry["properties"] = _mark(r)
            if r.status != "complete":
                failures += 1
                continue

            # Step 2 — TD-DFT. Extras/basis byte-match `oscpipe uvvis`, so a
            # prior uvvis run is a cache hit here and vice versa.
            method = args.method.lower()
            job, label = make_job(
                canonical,
                method,
                "6-31g**",
                0,
                1,
                "tddft",
                extras=f"n{args.nstates}",
                workflow_id=wf_id,
                notes=f"nstates={args.nstates}",
            )

            def _build_td(canonical=canonical, label=label, method=method):
                atoms = chem_smiles.embed_3d(canonical)
                return gaussian.write_com_tddft(
                    atoms,
                    method=method,
                    basis="6-31g**",
                    charge=0,
                    mult=1,
                    nstates=args.nstates,
                    nproc=settings.gaussian_nproc,
                    mem=settings.gaussian_mem,
                    label=label,
                    chk=f"{label}.chk",
                )

            r = runner.resolve(job, label, _build_td, wait=True)
            entry["tddft"] = _mark(r)
            if r.status != "complete":
                failures += 1
                continue

            # Step 3 — λ_h, the existing workflow, unmodified. It inserts its
            # own child workflow row; latest-by-smiles is safe single-process.
            rc = run_lambda_h(
                SimpleNamespace(smiles=canonical),
                settings,
                backend,
                conn,
                stdout=stdout,
                geometry_loader=geometry_loader,
            )
            child = conn.execute(
                "SELECT id, status FROM workflows "
                "WHERE kind='lambda_h' AND smiles=? ORDER BY id DESC LIMIT 1",
                (canonical,),
            ).fetchone()
            if child is not None:
                entry["lambda_h"] = {"workflow_id": child["id"], "status": child["status"]}
            if rc != 0:
                failures += 1
        except Exception as exc:
            # Continue-on-failure: one bad candidate (embed failure, SSH death
            # after the backend's reconnect-once gave up, …) must not kill a
            # weeks-long campaign.
            entry["error"] = str(exc)
            failures += 1
            print(f"screen: candidate failed: {exc}", file=stdout)

    db.update_workflow(
        conn,
        wf_id,
        "complete" if failures == 0 else "error",
        summary_json=json.dumps(
            {
                "total": len(smiles_list),
                "ok": len(smiles_list) - failures,
                "failed": failures,
                "candidates": progress,
            }
        ),
    )
    print(f"screen: {len(smiles_list) - failures}/{len(smiles_list)} ok", file=stdout)
    return 0 if failures == 0 else 1
