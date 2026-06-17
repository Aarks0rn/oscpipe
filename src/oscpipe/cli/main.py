"""oscpipe CLI entry point.

Subcommands:
    submit      Submit one molecule (Properties workflow). Fire-and-forget;
                polls once before exit so the LocalBackend's synchronous
                runs land in the DB inline.
    batch       Submit a list of SMILES (CSV with `smiles` column, or stdin).
    screen      Campaign: properties + uvvis + λ_h per candidate, sequential
                and blocking; resume by re-running (signature cache).
    fetch       Poll running jobs and pull completed logs back.
    status      Print the job table.
    export      Flat one-row-per-candidate CSV of stored results.
    reconcile   Poll backend for all pending/running jobs; mark losses.
    preflight   Workstation health check: g16, work_dir, DB.
    lambda      λ_reorg (4-point).
    uvvis       TDDFT excited states.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import sys
import time
from types import SimpleNamespace
from typing import IO

from ..chem import smiles as chem_smiles
from ..chem.geometry import geometry_hash, read_gaussian_log
from ..dft import gaussian
from ..run import JobRunner, _now, make_job
from ..run import db_dir as _db_dir
from ..settings import Settings, load
from ..store import db
from ..workflows.lambda_h import run_lambda_h
from ..workflows.oligomer_sweep import run_oligomer_sweep
from ..workflows.screen import run_screen as run_screen_workflow


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

    sc = sub.add_parser("screen", help="campaign: properties + uvvis + λ_h per candidate (blocking)")
    sc.add_argument("file", help='CSV with a "smiles" column, or "-" for stdin')
    sc.add_argument("--method", default="b3lyp")
    sc.add_argument("--basis", default="6-31g**")
    sc.add_argument("--nstates", type=int, default=10)

    lm = sub.add_parser("lambda", help="λ_reorg + Marcus rate (4-point)")
    lm.add_argument("smiles")

    ex = sub.add_parser("export", help="flat one-row-per-candidate CSV of stored results")
    ex.add_argument("--csv", required=True, help="output path")

    uv = sub.add_parser("uvvis", help="TDDFT excited states")
    uv.add_argument("smiles")
    uv.add_argument("--nstates", type=int, default=10)
    uv.add_argument("--method", default="b3lyp")
    uv.add_argument(
        "--from-log",
        dest="from_log",
        default=None,
        help=(
            "Reuse the final optimised geometry from this Gaussian .log "
            "(typically the neutral_opt.log produced by `oscpipe lambda`) "
            "instead of running a fresh RDKit 3D embed. Atom count and "
            "element multiset must match the SMILES."
        ),
    )

    og = sub.add_parser("oligomer", help="oligomer-length sweep + 1/n extrapolation")
    og.add_argument("repeat_unit", help="repeat-unit SMILES with two [*] attachment points")
    og.add_argument("--max-n", dest="max_n", type=int, default=3)
    og.add_argument("--nstates", type=int, default=10)

    st = sub.add_parser("status", help="print job table")
    st.add_argument(
        "--watch",
        type=int,
        default=0,
        metavar="N",
        help="refresh every N seconds (live); polls real PID liveness for running jobs",
    )
    st.add_argument("--running", action="store_true", help="show only pending/running jobs")
    sub.add_parser("preflight", help="workstation health check")
    sub.add_parser("reconcile", help="sync local DB with the backend")

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


def _submit_rc(r, stdout) -> int:
    """Map a fire-and-forget :class:`oscpipe.run.Resolved` to a CLI exit code.

    Shared tail of ``submit`` and ``uvvis`` (the cache-hit case is handled by
    each command, which prints its own message). ``complete`` → 0; an async
    backend that hasn't finished → 0 (user fetches later); a backend error → 1;
    a parse failure (log fetched but unparseable) → 2.
    """
    if r.status == "complete":
        return 0
    if r.status == "error":
        if r.log_path is None:
            print(f"job {r.job_id} ended with status=error", file=stdout)
            return 1
        return 2  # parse failure — _persist_by_kind already printed the reason
    return 0  # async backend; user fetches later


# ── submit ─────────────────────────────────────────────────────────────────


def run_submit(args, settings: Settings, backend, conn, *, stdout=None, workflow_id=None) -> int:
    if stdout is None:
        stdout = sys.stdout
    canonical, warnings = chem_smiles.canonicalise(args.smiles)
    for w in warnings:
        print(f"warning: {w}", file=stdout)

    job, label = make_job(
        canonical,
        args.method,
        args.basis,
        args.charge,
        args.mult,
        "properties",
        workflow_id=workflow_id,
    )

    def _build_com() -> str:
        atoms = chem_smiles.embed_3d(canonical)
        return gaussian.write_com_properties(
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

    runner = JobRunner(conn, settings, backend, stdout)
    r = runner.resolve(job, label, _build_com, wait=False)

    if r.cached:
        print(
            f"cache hit: job={r.job_id} HOMO={r.result.homo_ev:.3f} eV "
            f"LUMO={r.result.lumo_ev:.3f} eV gap={r.result.gap_ev:.3f} eV",
            file=stdout,
        )
        return 0
    return _submit_rc(r, stdout)


# ── uvvis ──────────────────────────────────────────────────────────────────


def run_uvvis(
    args,
    settings: Settings,
    backend,
    conn,
    *,
    stdout=None,
    geometry_loader=None,
) -> int:
    if stdout is None:
        stdout = sys.stdout
    if geometry_loader is None:
        geometry_loader = _atoms_from_log_validated
    canonical, warnings = chem_smiles.canonicalise(args.smiles)
    for w in warnings:
        print(f"warning: {w}", file=stdout)

    method = args.method.lower()
    from_log = getattr(args, "from_log", None)
    extras = f"n{args.nstates}"
    log_atoms = None
    if from_log:
        # Load eagerly so the geometry reaches the signature: a --from-log run
        # must never share a cache entry with a plain embed run of the same
        # SMILES (same molecule, different geometry).
        log_atoms = geometry_loader(from_log, canonical)
        print(f"uvvis: reusing geometry from {from_log}", file=stdout)
        extras += f"|geom{geometry_hash(log_atoms)}"
    job, label = make_job(
        canonical,
        method,
        "6-31g**",
        0,
        1,
        "tddft",
        extras=extras,
        notes=f"nstates={args.nstates}",
    )

    def _build_com() -> str:
        atoms = log_atoms if log_atoms is not None else chem_smiles.embed_3d(canonical)
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

    runner = JobRunner(conn, settings, backend, stdout)
    r = runner.resolve(job, label, _build_com, wait=False)

    if r.cached:
        print(f"cache hit: job={r.job_id} tddft (cached)", file=stdout)
        return 0
    return _submit_rc(r, stdout)


# ── lambda (4-point reorganisation energy) ─────────────────────────────────


def _atoms_from_log_validated(log_path: str, canonical_smiles: str):
    """Load the final optimised geometry from a Gaussian .log and sanity-check
    that its atom multiset matches what `canonical_smiles` would produce.

    Raises ``FileNotFoundError`` if the log is missing, ``ValueError`` if the
    log cannot be parsed by ASE or the atom multiset does not match the SMILES.
    Full connectivity is not checked — bond perception from coordinates is
    fragile for conjugated systems, so we settle for atom count + element
    multiset which catches the "wrong log file" failure mode cleanly.
    """
    from collections import Counter
    from pathlib import Path

    p = Path(log_path)
    if not p.exists():
        raise FileNotFoundError(f"--from-log: {log_path} does not exist")
    try:
        atoms = read_gaussian_log(log_path)
    except Exception as e:
        raise ValueError(f"--from-log: could not parse {log_path}: {e}") from e

    reference = chem_smiles.embed_3d(canonical_smiles)
    got = Counter(atoms.get_chemical_symbols())
    want = Counter(reference.get_chemical_symbols())
    if got != want:
        raise ValueError(
            f"--from-log: atom multiset mismatch between {log_path} and SMILES "
            f"{canonical_smiles!r}. log has {dict(got)}, SMILES wants {dict(want)}. "
            f"Is this the right .log for this molecule?"
        )
    return atoms


def run_lambda(
    args,
    settings: Settings,
    backend,
    conn,
    *,
    stdout=None,
    geometry_loader=None,
) -> int:
    """λ_h workflow — thin adapter over oscpipe.workflows.lambda_h.run_lambda_h."""
    return run_lambda_h(
        args, settings, backend, conn, stdout=stdout, geometry_loader=geometry_loader
    )


def run_oligomer(
    args,
    settings: Settings,
    backend,
    conn,
    *,
    stdout=None,
    geometry_loader=None,
) -> int:
    """Oligomer sweep — thin adapter over oscpipe.workflows.oligomer_sweep.run_oligomer_sweep."""
    return run_oligomer_sweep(
        args, settings, backend, conn, stdout=stdout, geometry_loader=geometry_loader
    )


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


# ── screen ─────────────────────────────────────────────────────────────────


def run_screen(args, settings: Settings, backend, conn, *, stdin=None, stdout=None) -> int:
    """Screen campaign — thin adapter over oscpipe.workflows.screen.run_screen."""
    if stdout is None:
        stdout = sys.stdout
    if stdin is None:
        stdin = sys.stdin
    smiles_list = _read_smiles(args.file, stdin)
    if not smiles_list:
        print("screen: no SMILES read", file=stdout)
        return 1
    return run_screen_workflow(smiles_list, args, settings, backend, conn, stdout=stdout)


# ── export ─────────────────────────────────────────────────────────────────


EXPORT_COLUMNS = [
    "smiles",
    "homo_ev",
    "lumo_ev",
    "gap_ev",
    "lambda_max_nm",
    "f_osc",
    "lambda_h_ev",
    "j_hole_ev",
    "marcus_rate_s1",
]


def run_export(args, settings: Settings, conn, *, stdout=None) -> int:
    """Flat one-row-per-candidate CSV of stored results (read-only, no backend)."""
    if stdout is None:
        stdout = sys.stdout
    smiles_list = db.list_candidate_smiles(conn)
    with open(args.csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=EXPORT_COLUMNS)
        writer.writeheader()
        for smi in smiles_list:
            row = db.candidate_summary(conn, smi)
            writer.writerow({k: ("" if row[k] is None else row[k]) for k in EXPORT_COLUMNS})
    print(f"export: {len(smiles_list)} candidates -> {args.csv}", file=stdout)
    return 0


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
    runner = JobRunner(conn, settings, backend, stdout)
    runner.rehydrate(rows)
    n_complete = n_error = n_still_running = 0
    for row in rows:
        status = backend.poll(row["ssh_jobid"])
        if status == "complete":
            final = runner.finish_completed(row)
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


_LIVE_LABELS = {"running": "ALIVE", "complete": "DONE", "error": "ERR", "unknown": "GONE"}


def _elapsed(started_at) -> str:
    """Human 'time since started_at' for the watch view; '—' if unparseable."""
    if not started_at:
        return "—"
    try:
        start = datetime.datetime.fromisoformat(started_at)
    except ValueError:
        return "—"
    secs = max(0, int((datetime.datetime.now() - start).total_seconds()))
    h, rem = divmod(secs, 3600)
    m, s = divmod(rem, 60)
    if h:
        return f"{h}h{m:02d}m"
    return f"{m}m{s:02d}s" if m else f"{s}s"


def _print_status_table(conn, backend, *, running_only: bool, stdout) -> None:
    """One snapshot of the job table. With ``backend``, running rows get a live
    column: poll the real PID (ALIVE/DONE/ERR/GONE) plus elapsed + workflow id."""
    where = "WHERE status IN ('pending', 'running')" if running_only else ""
    rows = list(
        conn.execute(
            "SELECT id, smiles, method, basis, charge, status, started_at, ssh_jobid, "
            f"workflow_id, homo_ev, gap_ev FROM v_jobs_with_results {where} ORDER BY id"
        )
    )
    if not rows:
        print("(no running jobs)" if running_only else "(no jobs)", file=stdout)
        return
    show_live = backend is not None
    # `status` is the raw job monitor (shows every job, incl. λ_h sub-jobs). The
    # charge column gives the orbital values context — a cation_opt / sp_cation row
    # is chg=1, so its HOMO/gap is not mistaken for a neutral result.
    live_hdr = f"{'live':<5}  {'elapsed':>8}  {'wf':>4}  " if show_live else ""
    print(
        f"{'id':>4}  {'status':<10}  {live_hdr}{'method/basis':<20}  "
        f"{'chg':>3}  {'HOMO':>8}  {'gap':>6}  smiles",
        file=stdout,
    )
    for r in rows:
        mb = f"{r['method']}/{r['basis']}"
        homo = f"{r['homo_ev']:.3f}" if r["homo_ev"] is not None else "—"
        gap = f"{r['gap_ev']:.3f}" if r["gap_ev"] is not None else "—"
        live_cols = ""
        if show_live:
            live = "—"
            if r["status"] == "running" and r["ssh_jobid"]:
                try:
                    live = _LIVE_LABELS.get(backend.poll(r["ssh_jobid"]), "?")
                except Exception:
                    live = "?"
            wf = r["workflow_id"] if r["workflow_id"] is not None else "—"
            live_cols = f"{live:<5}  {_elapsed(r['started_at']):>8}  {str(wf):>4}  "
        print(
            f"{r['id']:>4}  {r['status']:<10}  {live_cols}{mb:<20}  "
            f"{r['charge']:>3}  {homo:>8}  {gap:>6}  {r['smiles']}",
            file=stdout,
        )


def run_status(args, settings: Settings, conn, *, backend=None, stdout=None) -> int:
    if stdout is None:
        stdout = sys.stdout
    running_only = getattr(args, "running", False)
    # Live PID liveness needs the backend's jobid→log map populated for jobs this
    # process didn't submit; rehydrate from the running rows (no-op for LocalBackend).
    if backend is not None:
        inflight = list(conn.execute("SELECT * FROM jobs WHERE status = 'running'"))
        JobRunner(conn, settings, backend, stdout).rehydrate(inflight)

    watch = getattr(args, "watch", 0) or 0
    if watch <= 0:
        _print_status_table(conn, backend, running_only=running_only, stdout=stdout)
        return 0
    try:
        while True:
            stdout.write("\x1b[2J\x1b[H")  # clear screen + cursor home
            print(f"oscpipe status — refresh {watch}s — {_now()} (Ctrl-C to stop)", file=stdout)
            _print_status_table(conn, backend, running_only=running_only, stdout=stdout)
            stdout.flush()
            time.sleep(watch)
    except KeyboardInterrupt:
        return 0


# ── reconcile ──────────────────────────────────────────────────────────────


def run_reconcile(args, settings: Settings, backend, conn, *, stdout=None) -> int:
    if stdout is None:
        stdout = sys.stdout
    rows = list(conn.execute("SELECT * FROM jobs WHERE status IN ('pending','running')"))
    if not rows:
        print("reconcile: nothing to do", file=stdout)
        return 0
    runner = JobRunner(conn, settings, backend, stdout)
    runner.rehydrate(rows)
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
            runner.finish_completed(row)
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
        if args.cmd == "screen":
            return run_screen(args, settings, _make_backend(settings), conn)
        if args.cmd == "fetch":
            return run_fetch(args, settings, _make_backend(settings), conn)
        if args.cmd == "status":
            # Live watch polls real PID liveness, so it needs a backend; a one-shot
            # status stays DB-only (no SSH). Tolerate a misconfigured backend.
            backend = None
            if getattr(args, "watch", 0):
                try:
                    backend = _make_backend(settings)
                except ValueError:
                    backend = None
            return run_status(args, settings, conn, backend=backend)
        if args.cmd == "export":
            return run_export(args, settings, conn)
        if args.cmd == "reconcile":
            return run_reconcile(args, settings, _make_backend(settings), conn)
        if args.cmd == "preflight":
            return run_preflight(args, settings, _make_backend(settings), conn)
        if args.cmd == "uvvis":
            return run_uvvis(args, settings, _make_backend(settings), conn)
        if args.cmd == "lambda":
            return run_lambda(args, settings, _make_backend(settings), conn)
        if args.cmd == "oligomer":
            return run_oligomer(args, settings, _make_backend(settings), conn)
        raise NotImplementedError(f"cmd={args.cmd} not yet wired up")
    finally:
        conn.close()


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
