"""oscpipe CLI entry point.

Subcommands:
    submit      Submit one molecule (Properties workflow). Fire-and-forget;
                polls once before exit so the LocalBackend's synchronous
                runs land in the DB inline.
    batch       Submit a list of SMILES (CSV with `smiles` column, or stdin).
    fetch       Poll running jobs and pull completed logs back.
    status      Print the job table.
    reconcile   Poll backend for all pending/running jobs; mark losses.
    preflight   Workstation health check: g16, qstat, scratch_dir, DB.
    lambda      λ_reorg (4-point).
    uvvis       TDDFT excited states.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import sys
import tempfile
from collections.abc import Iterable
from pathlib import Path
from types import SimpleNamespace
from typing import IO

from ..chem import smiles as chem_smiles
from ..dft import gaussian
from ..settings import Settings, load
from ..store import db
from ..store.cache import signature


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="oscpipe")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("submit", help="submit one molecule")
    s.add_argument("smiles")
    s.add_argument("--method", default="b3lyp")
    s.add_argument("--basis", default="6-31g*")
    s.add_argument("--charge", type=int, default=0)
    s.add_argument("--mult", type=int, default=1)

    b = sub.add_parser("batch", help="submit a list of SMILES")
    b.add_argument("file", help='CSV with a "smiles" column, or "-" for stdin')
    b.add_argument("--method", default="b3lyp")
    b.add_argument("--basis", default="6-31g*")
    b.add_argument("--charge", type=int, default=0)
    b.add_argument("--mult", type=int, default=1)

    lm = sub.add_parser("lambda", help="λ_reorg + Marcus rate (4-point)")
    lm.add_argument("smiles")

    uv = sub.add_parser("uvvis", help="TDDFT excited states")
    uv.add_argument("smiles")
    uv.add_argument("--nstates", type=int, default=10)
    uv.add_argument("--method", default="b3lyp")

    sub.add_parser("status", help="print job table")
    sub.add_parser("preflight", help="workstation health check")
    sub.add_parser("reconcile", help="sync local DB with qstat")

    f = sub.add_parser("fetch", help="download completed logs + parse")
    f.add_argument("--job-id", type=int)

    return p


# ── helpers ────────────────────────────────────────────────────────────────


def _make_backend(settings: Settings):
    if settings.backend == "local":
        from ..dispatch.local import LocalBackend

        work = settings.local_work_dir or str(_db_dir(settings) / "work")
        return LocalBackend(work, exe=settings.gaussian_exe)
    if settings.backend == "ssh":
        from ..dispatch.ssh import SshBackend

        return SshBackend(settings)
    raise ValueError(f"unknown backend: {settings.backend!r}")


def _db_dir(settings: Settings) -> Path:
    return Path(settings.db_path).resolve().parent


def _log_dir(settings: Settings) -> Path:
    d = _db_dir(settings) / "logs"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _label(smiles_canonical: str, sig: str) -> str:
    slug = "".join(c if c.isalnum() else "_" for c in smiles_canonical)[:24]
    return f"{slug}_{sig}"


def _now() -> str:
    return datetime.datetime.now().isoformat(timespec="seconds")


def _remote_log_path(settings: Settings, label: str) -> str:
    base = (settings.remote_work_dir or "").rstrip("/")
    return f"{base}/{label}.log" if base else f"{label}.log"


def _rehydrate(backend, settings: Settings, rows: Iterable) -> None:
    """For SshBackend, repopulate jobid → remote_log_path from DB rows."""
    if not hasattr(backend, "_jobs"):
        return
    for r in rows:
        if not r["ssh_jobid"]:
            continue
        label = _label(r["smiles"], r["signature"])
        # Lambda workflow notes are plain identifiers used as label suffixes
        # (e.g. "neutral_opt", "cation_at_neutral_geom"). Metadata notes like
        # "nstates=10" contain "=" and must NOT be appended to the label.
        if r["notes"] and "=" not in r["notes"]:
            label = f"{label}_{r['notes']}"
        backend._jobs[r["ssh_jobid"]] = _remote_log_path(settings, label)


def _process_completion(conn, settings: Settings, backend, row, *, stdout: IO[str]) -> str:
    """Fetch log, parse, update DB for one row already known to be complete.

    Dispatches on row["job_kind"]: TDDFT rows parse excited states into
    spectra_json; everything else parses properties (HOMO/LUMO/gap/dipole/E).
    """
    label = _label(row["smiles"], row["signature"])
    if row["notes"] and "=" not in row["notes"]:
        label = f"{label}_{row['notes']}"
    local_log = backend.fetch_log(row["ssh_jobid"], label, str(_log_dir(settings)))
    try:
        if row["job_kind"] == "tddft":
            return _persist_tddft(conn, row, local_log, stdout)
        return _persist_properties(conn, row, local_log, stdout)
    except ValueError as exc:
        db.update_job_status(
            conn,
            row["id"],
            "error",
            completed_at=_now(),
            log_path=local_log,
            error_msg=str(exc),
        )
        print(f"job {row['id']} complete but parse failed: {exc}", file=stdout)
        return "error"


def _persist_properties(conn, row, local_log: str, stdout: IO[str]) -> str:
    props = gaussian.parse_properties(local_log)
    db.insert_result(
        conn,
        db.Result(
            job_id=row["id"],
            homo_ev=props.homo_ev,
            lumo_ev=props.lumo_ev,
            gap_ev=props.gap_ev,
            dipole_debye=props.dipole_debye,
            energy_ev=props.energy_ev,
        ),
    )
    db.update_job_status(conn, row["id"], "complete", completed_at=_now(), log_path=local_log)
    print(
        f"complete: job={row['id']} HOMO={props.homo_ev:.3f} eV "
        f"LUMO={props.lumo_ev:.3f} eV gap={props.gap_ev:.3f} eV",
        file=stdout,
    )
    return "complete"


def _persist_tddft(conn, row, local_log: str, stdout: IO[str]) -> str:
    import json

    states = gaussian.parse_excited_states(local_log)
    if not states:
        raise ValueError(f"{local_log}: no excited states parsed")
    spectra = [
        {
            "n": s.n,
            "energy_ev": s.energy_ev,
            "wavelength_nm": s.wavelength_nm,
            "f": s.oscillator_strength,
        }
        for s in states
    ]
    db.insert_result(
        conn,
        db.Result(job_id=row["id"], spectra_json=json.dumps(spectra)),
    )
    db.update_job_status(conn, row["id"], "complete", completed_at=_now(), log_path=local_log)
    bright = max(states, key=lambda s: s.oscillator_strength)
    print(
        f"complete: job={row['id']} tddft n={len(states)} "
        f"brightest λ={bright.wavelength_nm:.1f} nm f={bright.oscillator_strength:.3f}",
        file=stdout,
    )
    return "complete"


# ── submit ─────────────────────────────────────────────────────────────────


def _submit_job_and_poll_once(
    conn,
    backend,
    *,
    job: db.Job,
    com_text: str,
    label: str,
    stdout: IO[str],
) -> tuple[int, str]:
    """Insert pending row, hand .com to backend, persist running state, poll once.

    Returns (job_id, status_after_one_poll). Synchronous backends report
    'complete' here; async backends typically return 'running'.
    """
    job_id = db.insert_job(conn, job)
    with tempfile.TemporaryDirectory() as tmp:
        com_path = Path(tmp) / f"{label}.com"
        com_path.write_text(com_text)
        ssh_jobid = backend.submit(str(com_path), label)
    db.update_job_status(conn, job_id, "running", started_at=_now(), ssh_jobid=ssh_jobid)
    print(f"submitted: job={job_id} jobid={ssh_jobid} label={label}", file=stdout)
    return job_id, backend.poll(ssh_jobid)


def run_submit(args, settings: Settings, backend, conn, *, stdout=None, workflow_id=None) -> int:
    if stdout is None:
        stdout = sys.stdout
    canonical, warnings = chem_smiles.canonicalise(args.smiles)
    for w in warnings:
        print(f"warning: {w}", file=stdout)

    sig = signature(canonical, args.method, args.basis, args.charge, args.mult)

    hit = db.find_complete_by_signature(conn, sig)
    if hit is not None:
        print(
            f"cache hit: job={hit['id']} HOMO={hit['homo_ev']:.3f} eV "
            f"LUMO={hit['lumo_ev']:.3f} eV gap={hit['gap_ev']:.3f} eV",
            file=stdout,
        )
        return 0

    atoms = chem_smiles.embed_3d(canonical)
    label = _label(canonical, sig)
    com_text = gaussian.write_com_properties(
        atoms,
        method=args.method,
        basis=args.basis,
        charge=args.charge,
        mult=args.mult,
        nproc=settings.gaussian_nproc,
        mem=settings.gaussian_mem,
        label=label,
        chk=f"{label}.chk",
    )

    job_id, status = _submit_job_and_poll_once(
        conn,
        backend,
        job=db.Job(
            id=None,
            signature=sig,
            smiles=canonical,
            method=args.method,
            basis=args.basis,
            charge=args.charge,
            mult=args.mult,
            job_kind="properties",
            status="pending",
            submitted_at=_now(),
            workflow_id=workflow_id,
        ),
        com_text=com_text,
        label=label,
        stdout=stdout,
    )

    if status not in ("complete", "error"):
        return 0  # async backend; user fetches later

    row = conn.execute("SELECT * FROM jobs WHERE id = ?", (job_id,)).fetchone()
    if status == "error":
        db.update_job_status(
            conn, job_id, "error", completed_at=_now(), error_msg="backend reported error"
        )
        print(f"job {job_id} ended with status=error", file=stdout)
        return 1
    final = _process_completion(conn, settings, backend, row, stdout=stdout)
    return 0 if final == "complete" else 2


# ── uvvis ──────────────────────────────────────────────────────────────────


def run_uvvis(args, settings: Settings, backend, conn, *, stdout=None) -> int:
    if stdout is None:
        stdout = sys.stdout
    canonical, warnings = chem_smiles.canonicalise(args.smiles)
    for w in warnings:
        print(f"warning: {w}", file=stdout)

    method = args.method.lower()
    sig = signature(
        canonical,
        method,
        "6-31g**",
        0,
        1,
        job_kind="tddft",
        extras=f"n{args.nstates}",
    )

    hit = db.find_complete_by_signature(conn, sig)
    if hit is not None:
        spectra = conn.execute(
            "SELECT spectra_json FROM results WHERE job_id = ?", (hit["id"],)
        ).fetchone()
        print(
            f"cache hit: job={hit['id']} tddft spectra={spectra['spectra_json'][:60]}...",
            file=stdout,
        )
        return 0

    atoms = chem_smiles.embed_3d(canonical)
    label = _label(canonical, sig)
    com_text = gaussian.write_com_tddft(
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

    job_id, status = _submit_job_and_poll_once(
        conn,
        backend,
        job=db.Job(
            id=None,
            signature=sig,
            smiles=canonical,
            method=method,
            basis="6-31g**",
            charge=0,
            mult=1,
            job_kind="tddft",
            status="pending",
            submitted_at=_now(),
            notes=f"nstates={args.nstates}",
        ),
        com_text=com_text,
        label=label,
        stdout=stdout,
    )

    if status not in ("complete", "error"):
        return 0

    row = conn.execute("SELECT * FROM jobs WHERE id = ?", (job_id,)).fetchone()
    if status == "error":
        db.update_job_status(
            conn, job_id, "error", completed_at=_now(), error_msg="backend reported error"
        )
        print(f"job {job_id} ended with status=error", file=stdout)
        return 1
    final = _process_completion(conn, settings, backend, row, stdout=stdout)
    return 0 if final == "complete" else 2


# ── lambda (4-point reorganisation energy) ─────────────────────────────────


def _build_pi_stack_dimer(atoms, stack_distance: float = 3.5):
    """Return ASE Atoms: two copies of atoms face-to-face, second shifted +z by stack_distance Å."""
    second = atoms.copy()
    second.positions[:, 2] += stack_distance
    return atoms + second


def _wait_for_job(backend, ssh_jobid: str, interval: int) -> str:
    """Block until backend reports a terminal state for the job.

    Ctrl-C detaches cleanly — job stays running on the workstation.
    Run `oscpipe fetch` later to pull results.
    """
    import time

    try:
        while True:
            s = backend.poll(ssh_jobid)
            if s in ("complete", "error", "unknown"):
                return s
            time.sleep(interval)
    except KeyboardInterrupt:
        return "detached"


def _submit_and_wait(
    conn,
    settings: Settings,
    backend,
    *,
    job: db.Job,
    com_text: str,
    label: str,
    stdout: IO[str],
) -> tuple[int, str, str | None]:
    """Submit a job, block until terminal, fetch log. Returns (job_id, status, log_path)."""
    job_id = db.insert_job(conn, job)
    with tempfile.TemporaryDirectory() as tmp:
        com_path = Path(tmp) / f"{label}.com"
        com_path.write_text(com_text)
        ssh_jobid = backend.submit(str(com_path), label)
    db.update_job_status(conn, job_id, "running", started_at=_now(), ssh_jobid=ssh_jobid)
    print(f"submitted: job={job_id} jobid={ssh_jobid} label={label}", file=stdout)

    status = _wait_for_job(backend, ssh_jobid, settings.poll_interval_seconds)
    if status == "detached":
        # Job still running on workstation — leave status=running in DB.
        print(
            f"detached: job={job_id} still running — run `oscpipe fetch` to pull results",
            file=stdout,
        )
        return job_id, status, None
    if status != "complete":
        db.update_job_status(
            conn,
            job_id,
            "error",
            completed_at=_now(),
            error_msg=f"backend returned {status}",
        )
        return job_id, status, None
    log_path = backend.fetch_log(ssh_jobid, label, str(_log_dir(settings)))
    db.update_job_status(conn, job_id, "complete", completed_at=_now(), log_path=log_path)
    return job_id, status, log_path


def _ase_read_gaussian(log_path: str):
    import ase.io

    return ase.io.read(log_path, format="gaussian-out")


def run_lambda(
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
        geometry_loader = _ase_read_gaussian

    import json

    canonical, warnings = chem_smiles.canonicalise(args.smiles)
    for w in warnings:
        print(f"warning: {w}", file=stdout)

    initial = chem_smiles.embed_3d(canonical)
    method, basis = "b3lyp", "6-31g**"

    wf_id = db.insert_workflow(conn, "lambda_h", canonical, _now())
    print(f"workflow: id={wf_id} kind=lambda_h smiles={canonical}", file=stdout)

    def _opt_job(charge, mult, suffix):
        sig = signature(canonical, method, basis, charge, mult)
        label = f"{_label(canonical, sig)}_{suffix}"
        com_text = gaussian.write_com_properties(
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
        return _submit_and_wait(
            conn, settings, backend, job=job, com_text=com_text, label=label, stdout=stdout
        )

    # 1+2: optimise neutral and cation.
    _, n_status, n_log = _opt_job(0, 1, "neutral_opt")
    _, c_status, c_log = _opt_job(1, 2, "cation_opt")
    if n_status != "complete" or c_status != "complete":
        db.update_workflow(
            conn,
            wf_id,
            "error",
            summary_json=json.dumps({"stage": "opt", "neutral": n_status, "cation": c_status}),
        )
        print("lambda: opt stage failed", file=stdout)
        return 1
    e_neut_opt = gaussian.parse_properties(n_log).energy_ev
    e_cat_opt = gaussian.parse_properties(c_log).energy_ev
    neutral_geom = geometry_loader(n_log)
    cation_geom = geometry_loader(c_log)

    # 3+4: single points at the *other* geometry.
    def _sp_job(atoms, charge, mult, kind, suffix):
        sig = signature(canonical, method, basis, charge, mult, extras=suffix)
        label = f"{_label(canonical, sig)}_{suffix}"
        com_text = gaussian.write_com_sp(
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
        return _submit_and_wait(
            conn, settings, backend, job=job, com_text=com_text, label=label, stdout=stdout
        )

    _, n_at_c_status, n_at_c_log = _sp_job(
        cation_geom, 0, 1, "sp_neutral", "neutral_at_cation_geom"
    )
    _, c_at_n_status, c_at_n_log = _sp_job(
        neutral_geom, 1, 2, "sp_cation", "cation_at_neutral_geom"
    )
    if n_at_c_status != "complete" or c_at_n_status != "complete":
        db.update_workflow(
            conn,
            wf_id,
            "error",
            summary_json=json.dumps(
                {"stage": "sp", "n_at_c": n_at_c_status, "c_at_n": c_at_n_status}
            ),
        )
        print("lambda: sp stage failed", file=stdout)
        return 1

    e_neut_at_cat = gaussian.parse_properties(n_at_c_log).energy_ev
    e_cat_at_neut = gaussian.parse_properties(c_at_n_log).energy_ev

    lambda_h = (e_cat_at_neut - e_cat_opt) + (e_neut_at_cat - e_neut_opt)
    print(f"lambda_h = {lambda_h:.4f} eV", file=stdout)

    # ── job 5: dimer SP → J_hole → Marcus rate ──────────────────────────────
    from ..analysis.marcus import marcus_rate as _marcus_rate

    dimer_atoms = _build_pi_stack_dimer(neutral_geom)
    sig_dimer = signature(canonical, method, basis, 0, 1, extras="dimer_sp")
    label_dimer = f"{_label(canonical, sig_dimer)}_dimer_sp"
    com_dimer = gaussian.write_com_sp(
        dimer_atoms,
        method=method,
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
        method=method,
        basis=basis,
        charge=0,
        mult=1,
        job_kind="sp_dimer",
        status="pending",
        submitted_at=_now(),
        workflow_id=wf_id,
        notes="dimer_sp",
    )
    _, dimer_status, dimer_log = _submit_and_wait(
        conn, settings, backend, job=job_dimer, com_text=com_dimer, label=label_dimer, stdout=stdout
    )

    j_hole_ev: float | None = None
    marcus: float | None = None
    if dimer_status == "complete" and dimer_log:
        homo, hm1, _lumo, _lp1 = gaussian.parse_dimer_orbitals(dimer_log)
        j_hole_ev = abs(homo - hm1) / 2.0
        marcus = _marcus_rate(lambda_h, j_hole_ev, delta_g_ev=0.0)
        print(f"J_hole = {j_hole_ev:.4f} eV", file=stdout)
        print(f"marcus_rate = {marcus:.4e} s^-1", file=stdout)
    else:
        print(f"dimer SP {dimer_status} — J/marcus not available", file=stdout)

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


# ── preflight ──────────────────────────────────────────────────────────────


def run_preflight(args, settings: Settings, backend, conn, *, stdout=None) -> int:
    if stdout is None:
        stdout = sys.stdout

    try:
        n = conn.execute("SELECT count(*) FROM jobs").fetchone()[0]
        print(f"db: ok ({n} jobs)", file=stdout)
    except Exception as exc:
        print(f"db: ERROR — {exc}", file=stdout)
        return 1

    if not hasattr(backend, "preflight"):
        print("backend: no preflight checks (local dev mode)", file=stdout)
        return 0

    checks = backend.preflight()
    all_ok = True
    for name, ok, msg in checks:
        marker = "ok" if ok else "FAIL"
        print(f"{name}: [{marker}] {msg}", file=stdout)
        if not ok:
            all_ok = False

    return 0 if all_ok else 1


# ── batch ──────────────────────────────────────────────────────────────────


def _read_smiles(path: str, stdin: IO[str]) -> list[str]:
    """CSV with header 'smiles' (and optional method/basis/charge/mult), or
    a plain one-SMILES-per-line stream when `path` == '-' or has no header.
    """
    if path == "-":
        text = stdin.read()
        return [ln.strip() for ln in text.splitlines() if ln.strip()]
    with open(path) as f:
        sample = f.read(2048)
        f.seek(0)
        first = (sample.lower().splitlines() or [""])[0]
        if "smiles" in first:
            reader = csv.DictReader(f)
            return [row["smiles"].strip() for row in reader if row.get("smiles", "").strip()]
        return [ln.strip() for ln in f if ln.strip()]


def run_batch(args, settings: Settings, backend, conn, *, stdin=None, stdout=None) -> int:
    if stdout is None:
        stdout = sys.stdout
    if stdin is None:
        stdin = sys.stdin
    smiles_list = _read_smiles(args.file, stdin)
    if not smiles_list:
        print("batch: no SMILES read", file=stdout)
        return 1

    wf_id = db.insert_workflow(conn, "batch", f"batch:{len(smiles_list)}", _now(), status="running")
    print(f"batch: workflow_id={wf_id} submitting {len(smiles_list)} molecules", file=stdout)

    failures = 0
    for smi in smiles_list:
        one = SimpleNamespace(
            smiles=smi,
            method=args.method,
            basis=args.basis,
            charge=args.charge,
            mult=args.mult,
        )
        rc = run_submit(one, settings, backend, conn, stdout=stdout, workflow_id=wf_id)
        if rc != 0:
            failures += 1

    final_status = "complete" if failures == 0 else "error"
    import json

    db.update_workflow(
        conn,
        wf_id,
        final_status,
        summary_json=json.dumps(
            {
                "total": len(smiles_list),
                "ok": len(smiles_list) - failures,
                "failed": failures,
            }
        ),
    )
    print(f"batch: {len(smiles_list) - failures}/{len(smiles_list)} ok", file=stdout)
    return 0 if failures == 0 else 1


# ── fetch ──────────────────────────────────────────────────────────────────


def run_fetch(args, settings: Settings, backend, conn, *, stdout=None) -> int:
    if stdout is None:
        stdout = sys.stdout
    if args.job_id is not None:
        rows = list(
            conn.execute(
                "SELECT * FROM jobs WHERE id = ? AND status IN ('pending','running')",
                (args.job_id,),
            )
        )
    else:
        rows = list(
            conn.execute(
                "SELECT * FROM jobs WHERE status IN ('pending','running') AND ssh_jobid IS NOT NULL"
            )
        )
    if not rows:
        print("fetch: no running jobs", file=stdout)
        return 0
    _rehydrate(backend, settings, rows)
    n_complete = n_error = n_still_running = 0
    for row in rows:
        status = backend.poll(row["ssh_jobid"])
        if status == "complete":
            final = _process_completion(conn, settings, backend, row, stdout=stdout)
            if final == "complete":
                n_complete += 1
            else:
                n_error += 1
        elif status == "error":
            db.update_job_status(
                conn,
                row["id"],
                "error",
                completed_at=_now(),
                error_msg="backend reported error",
            )
            print(f"job {row['id']}: error", file=stdout)
            n_error += 1
        else:
            n_still_running += 1
    print(
        f"fetch: {n_complete} complete, {n_error} error, {n_still_running} still running",
        file=stdout,
    )
    return 0


# ── status ─────────────────────────────────────────────────────────────────


def run_status(args, settings: Settings, conn, *, stdout=None) -> int:
    if stdout is None:
        stdout = sys.stdout
    rows = list(
        conn.execute(
            "SELECT id, smiles, method, basis, status, submitted_at, "
            "homo_ev, gap_ev FROM v_jobs_with_results ORDER BY id"
        )
    )
    if not rows:
        print("(no jobs)", file=stdout)
        return 0
    header = f"{'id':>4}  {'status':<10}  {'method/basis':<20}  {'HOMO':>8}  {'gap':>6}  smiles"
    print(header, file=stdout)
    for r in rows:
        mb = f"{r['method']}/{r['basis']}"
        homo = f"{r['homo_ev']:.3f}" if r["homo_ev"] is not None else "—"
        gap = f"{r['gap_ev']:.3f}" if r["gap_ev"] is not None else "—"
        print(
            f"{r['id']:>4}  {r['status']:<10}  {mb:<20}  {homo:>8}  {gap:>6}  {r['smiles']}",
            file=stdout,
        )
    return 0


# ── reconcile ──────────────────────────────────────────────────────────────


def run_reconcile(args, settings: Settings, backend, conn, *, stdout=None) -> int:
    if stdout is None:
        stdout = sys.stdout
    rows = list(conn.execute("SELECT * FROM jobs WHERE status IN ('pending','running')"))
    if not rows:
        print("reconcile: nothing to do", file=stdout)
        return 0
    _rehydrate(backend, settings, rows)
    promoted = lost = 0
    for row in rows:
        if not row["ssh_jobid"]:
            db.update_job_status(
                conn,
                row["id"],
                "error",
                completed_at=_now(),
                error_msg="lost: no ssh_jobid recorded",
            )
            lost += 1
            continue
        status = backend.poll(row["ssh_jobid"])
        if status == "complete":
            _process_completion(conn, settings, backend, row, stdout=stdout)
            promoted += 1
        elif status == "error":
            db.update_job_status(
                conn,
                row["id"],
                "error",
                completed_at=_now(),
                error_msg="backend reported error",
            )
            promoted += 1
        elif status == "unknown":
            db.update_job_status(
                conn,
                row["id"],
                "error",
                completed_at=_now(),
                error_msg="lost: backend has no record of this jobid",
            )
            lost += 1
    print(f"reconcile: {promoted} promoted, {lost} marked lost", file=stdout)
    return 0


# ── entry ──────────────────────────────────────────────────────────────────


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    settings = load()
    conn = db.open(settings.db_path)
    try:
        if args.cmd == "submit":
            return run_submit(args, settings, _make_backend(settings), conn)
        if args.cmd == "batch":
            return run_batch(args, settings, _make_backend(settings), conn)
        if args.cmd == "fetch":
            return run_fetch(args, settings, _make_backend(settings), conn)
        if args.cmd == "status":
            return run_status(args, settings, conn)
        if args.cmd == "reconcile":
            return run_reconcile(args, settings, _make_backend(settings), conn)
        if args.cmd == "preflight":
            return run_preflight(args, settings, _make_backend(settings), conn)
        if args.cmd == "uvvis":
            return run_uvvis(args, settings, _make_backend(settings), conn)
        if args.cmd == "lambda":
            return run_lambda(args, settings, _make_backend(settings), conn)
        raise NotImplementedError(f"cmd={args.cmd} not yet wired up")
    finally:
        conn.close()


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
