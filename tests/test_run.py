"""Tests for the JobRunner.resolve seam (cache → submit → persist-by-kind).

These exercise behaviour introduced when the single-shot CLI path and the λ_h
workflow were unified onto `JobRunner.resolve`:
  - a signature cache hit short-circuits without rebuilding the .com,
  - `need_log` recomputes when a cached job's local .log is gone,
  - persist-by-kind is total: sp_neutral/sp_cation get a results row, sp_dimer
    does not (its J_hole lives in the workflow summary),
  - re-running `oscpipe lambda` reuses completed sub-jobs from the cache.
"""

from __future__ import annotations

import argparse
import io
import json
from pathlib import Path

from _helpers import StubBackend, synth_log

from oscpipe.cli.main import run_fetch, run_lambda
from oscpipe.dft.gaussian import HARTREE_TO_EV
from oscpipe.run import JobRunner, JobSpec
from oscpipe.settings import Settings
from oscpipe.store import db
from oscpipe.store.cache import signature


def _settings(tmp_path):
    return Settings(
        backend="local",
        db_path=str(tmp_path / "results.db"),
        gaussian_nproc=2,
        gaussian_mem="2GB",
        poll_interval_seconds=0,
    )


def _job(sig, **overrides):
    base = dict(
        id=None,
        signature=sig,
        smiles="c1ccccc1",
        method="b3lyp",
        basis="6-31g**",
        charge=0,
        mult=1,
        job_kind="properties",
        status="pending",
        submitted_at="2026-05-29T00:00:00",
    )
    base.update(overrides)
    return db.Job(**base)


# λ_h per-suffix energies (match test_cli_lambda so λ_h ≠ 0, marcus stays finite).
_LAMBDA_ENERGIES_HA = {
    "neutral_opt": -154.0,
    "cation_opt": -153.5,
    "neutral_at_cation_geom": -153.95,
    "cation_at_neutral_geom": -153.45,
}


def _log_for_lambda(label: str) -> str:
    for suffix, energy in _LAMBDA_ENERGIES_HA.items():
        if label.endswith(suffix):
            return synth_log(energy_ha=energy)
    return synth_log()


# ── resolve: cache hit ───────────────────────────────────────────────────────


def test_resolve_cache_hit_returns_result_without_rebuilding(tmp_path):
    s = _settings(tmp_path)
    conn = db.open(s.db_path)
    sig = signature("c1ccccc1", "b3lyp", "6-31g**", 0, 1)
    jid = db.insert_job(conn, _job(sig, status="complete", log_path="/seeded.log"))
    db.insert_result(
        conn, db.Result(job_id=jid, homo_ev=-6.0, lumo_ev=-1.0, gap_ev=5.0, energy_ev=-4000.0)
    )

    runner = JobRunner(conn, s, StubBackend(), io.StringIO())
    built: list[int] = []

    def _build():
        built.append(1)
        return "com"

    r = runner.resolve(_job(sig), "lbl", _build, wait=False)

    assert r.cached is True
    assert r.result.energy_ev == -4000.0
    assert built == []  # .com was not rebuilt
    assert conn.execute("SELECT count(*) FROM jobs").fetchone()[0] == 1  # no new submit


def test_resolve_need_log_recomputes_when_cached_log_missing(tmp_path):
    s = _settings(tmp_path)
    conn = db.open(s.db_path)
    sig = signature("c1ccccc1", "b3lyp", "6-31g**", 0, 1)
    jid = db.insert_job(conn, _job(sig, status="complete", log_path=str(tmp_path / "gone.log")))
    db.insert_result(
        conn, db.Result(job_id=jid, homo_ev=-6.0, lumo_ev=-1.0, gap_ev=5.0, energy_ev=-4000.0)
    )

    backend = StubBackend(log_text=synth_log(energy_ha=-100.0))
    backend.poll_return = "complete"
    runner = JobRunner(conn, s, backend, io.StringIO())
    built: list[int] = []

    def _build():
        built.append(1)
        return "com"

    r = runner.resolve(_job(sig), "lbl", _build, wait=True, need_log=True)

    assert built == [1]  # cached .log gone → recomputed
    assert r.cached is False
    assert r.status == "complete"
    assert conn.execute("SELECT count(*) FROM jobs").fetchone()[0] == 2


# ── persist-by-kind totality (through λ_h) ───────────────────────────────────


def test_lambda_persists_sp_results_but_not_dimer(tmp_path):
    s = _settings(tmp_path)
    backend = StubBackend(log_provider=_log_for_lambda)
    backend.poll_return = "complete"
    conn = db.open(s.db_path)
    from oscpipe.chem.smiles import embed_3d

    initial = embed_3d("c1ccccc1")
    rc = run_lambda(
        argparse.Namespace(smiles="c1ccccc1"), s, backend, conn, geometry_loader=lambda _: initial
    )
    assert rc == 0

    for kind in ("sp_neutral", "sp_cation"):
        jrow = conn.execute("SELECT id FROM jobs WHERE job_kind = ?", (kind,)).fetchone()
        res = conn.execute(
            "SELECT energy_ev FROM results WHERE job_id = ?", (jrow["id"],)
        ).fetchone()
        assert res is not None, f"{kind} should have a results row"
        assert res["energy_ev"] is not None

    drow = conn.execute("SELECT id FROM jobs WHERE job_kind = 'sp_dimer'").fetchone()
    assert conn.execute("SELECT 1 FROM results WHERE job_id = ?", (drow["id"],)).fetchone() is None


def test_lambda_rerun_reuses_cache_no_new_jobs(tmp_path):
    s = _settings(tmp_path)
    backend = StubBackend(log_provider=_log_for_lambda)
    backend.poll_return = "complete"
    conn = db.open(s.db_path)
    from oscpipe.chem.smiles import embed_3d

    initial = embed_3d("c1ccccc1")
    args = argparse.Namespace(smiles="c1ccccc1")

    assert run_lambda(args, s, backend, conn, geometry_loader=lambda _: initial) == 0
    n_after_first = conn.execute("SELECT count(*) FROM jobs").fetchone()[0]
    assert n_after_first == 5

    # Second run: every sub-job signature hits the cache → no new jobs submitted.
    assert run_lambda(args, s, backend, conn, geometry_loader=lambda _: initial) == 0
    assert conn.execute("SELECT count(*) FROM jobs").fetchone()[0] == 5
    # Two workflows recorded, both complete.
    wfs = list(
        conn.execute(
            "SELECT status, summary_json FROM workflows WHERE kind = 'lambda_h' ORDER BY id"
        )
    )
    assert len(wfs) == 2
    assert all(w["status"] == "complete" for w in wfs)
    # The cached re-run reproduces the same λ_h, not just the same job count.
    second = json.loads(wfs[1]["summary_json"])
    assert abs(second["lambda_h_ev"] - 2.72114) < 1e-2
    assert second["j_hole_ev"] is not None


# ── persist-by-kind on the resume (fetch) path ───────────────────────────────


def test_fetch_resumes_dimer_without_bogus_results_row(tmp_path):
    s = _settings(tmp_path)
    conn = db.open(s.db_path)
    sig = signature("c1ccccc1", "wb97xd", "6-31g**", 0, 1, extras="dimer_sp")
    jid = db.insert_job(
        conn, _job(sig, method="wb97xd", job_kind="sp_dimer", status="running", ssh_jobid="d-1")
    )

    backend = StubBackend()  # default SYNTH_LOG has orbitals — would mis-persist as properties
    backend.poll_return = "complete"
    run_fetch(argparse.Namespace(job_id=None), s, backend, conn)

    row = conn.execute("SELECT status FROM jobs WHERE id = ?", (jid,)).fetchone()
    assert row["status"] == "complete"
    assert conn.execute("SELECT 1 FROM results WHERE job_id = ?", (jid,)).fetchone() is None


# ── Properties browser filter (mirrors app/pages/1_Properties.py) ─────────────


def test_properties_browser_hides_lambda_subjobs(tmp_path):
    """Total-persist gives λ_h opt jobs (job_kind='properties') a results row; the
    Properties browser must hide those workflow internals while still showing
    standalone and batch property jobs."""
    s = _settings(tmp_path)
    conn = db.open(s.db_path)
    lam_wf = db.insert_workflow(conn, "lambda_h", "c1ccccc1", "2026-05-29T00:00:00")
    batch_wf = db.insert_workflow(conn, "batch", "batch:1", "2026-05-29T00:00:00")

    def _complete_prop(smiles, sig, workflow_id):
        jid = db.insert_job(
            conn, _job(sig, smiles=smiles, status="complete", workflow_id=workflow_id)
        )
        db.insert_result(
            conn,
            db.Result(job_id=jid, homo_ev=-6.0, lumo_ev=-1.0, gap_ev=5.0, energy_ev=-4000.0),
        )
        return jid

    lam_job = _complete_prop("c1ccccc1", "sig-lam", lam_wf)
    standalone = _complete_prop("c1ccsc1", "sig-std", None)
    batch_job = _complete_prop("c1ccoc1", "sig-bat", batch_wf)

    rows = conn.execute(
        "SELECT id FROM v_jobs_with_results WHERE job_kind = 'properties' "
        "AND (workflow_id IS NULL OR workflow_id NOT IN "
        "(SELECT id FROM workflows WHERE kind = 'lambda_h')) "
        "ORDER BY id"
    ).fetchall()
    ids = {r["id"] for r in rows}
    assert standalone in ids  # standalone submit shown
    assert batch_job in ids  # batch job (has workflow_id, but kind=batch) shown
    assert lam_job not in ids  # λ_h sub-job hidden from the properties browser


def test_make_job_matches_existing_signature_convention():
    """make_job must reproduce the exact pre-refactor signature/label forms —
    a drift here would orphan every live cache entry."""
    from oscpipe.run import label, make_job
    from oscpipe.store.cache import signature

    canon = "c1ccccc1"
    job, lbl = make_job(canon, "b3lyp", "6-31g*", 0, 1, "properties")
    assert job.signature == signature(canon, "b3lyp", "6-31g*", 0, 1)
    assert lbl == label(canon, job.signature)
    assert job.notes is None

    job2, _ = make_job(
        canon, "b3lyp", "6-31g**", 0, 1, "tddft", extras="n10", notes="nstates=10"
    )
    assert job2.signature == signature(
        canon, "b3lyp", "6-31g**", 0, 1, job_kind="tddft", extras="n10"
    )
    assert job2.notes == "nstates=10"


# ── resolve: reattach to an in-flight job instead of resubmitting ────────────


def test_resolve_reattaches_to_inflight_job_no_resubmit(tmp_path):
    """A prior detached run left a 'running' row with this signature + ssh_jobid.
    Re-running must poll that job, not qsub a duplicate (orphan prevention)."""
    s = _settings(tmp_path)
    conn = db.open(s.db_path)
    sig = signature("c1ccccc1", "b3lyp", "6-31g**", 0, 1)
    existing_id = db.insert_job(
        conn, _job(sig, status="running", ssh_jobid="71229", notes="n2_opt")
    )

    backend = StubBackend(log_text=synth_log())
    backend.poll_return = "complete"
    runner = JobRunner(conn, s, backend, io.StringIO())

    def _build():
        raise AssertionError("build_com must not run on reattach")

    r = runner.resolve(_job(sig, notes="n2_opt"), "lbl", _build, wait=True)

    assert r.job_id == existing_id          # reused the existing row
    assert ("submit", "lbl") not in backend.calls and not any(
        c[0] == "submit" for c in backend.calls
    )                                       # never qsub'd a duplicate
    assert conn.execute("SELECT count(*) FROM jobs").fetchone()[0] == 1
    assert r.status == "complete"           # polled the reattached job to completion
    conn.close()


def test_resolve_submits_when_only_terminal_rows_exist(tmp_path):
    """An 'error' row is terminal — it must NOT be reattached (resubmit instead)."""
    s = _settings(tmp_path)
    conn = db.open(s.db_path)
    sig = signature("c1ccccc1", "b3lyp", "6-31g**", 0, 1)
    db.insert_job(conn, _job(sig, status="error", ssh_jobid="71231", notes="n3_opt"))

    backend = StubBackend(log_text=synth_log())
    backend.poll_return = "complete"
    runner = JobRunner(conn, s, backend, io.StringIO())

    built = []

    def _build():
        built.append(True)
        return "stub com"

    runner.resolve(_job(sig, notes="n3_opt"), "lbl", _build, wait=True)

    assert built == [True]                  # built + submitted a fresh job
    assert any(c[0] == "submit" for c in backend.calls)
    assert conn.execute("SELECT count(*) FROM jobs").fetchone()[0] == 2
    conn.close()


# ── preopt accepts the 'steps exceeded' termination as a usable seed ──────────

_STEPS_EXCEEDED_LOG = (
    " GradGradGrad\n"
    "    -- Number of steps exceeded,  NStep= 300\n"
    " Error termination via Lnk1e in /usr/local/g16/l9999.exe\n"
)


def test_resolve_preopt_accepts_steps_exceeded_as_seed(tmp_path):
    # A PM6 pre-opt is a geometry SEED: floppy oligomers oscillate just above even
    # the loose target and g16 Error-terminates 'Number of steps exceeded', but the
    # last geometry is a fine seed for the B3LYP refine — so it must resolve to
    # 'complete', not 'error' (job 71251).
    s = _settings(tmp_path)
    conn = db.open(s.db_path)
    sig = signature("c1ccccc1", "pm6", "", 0, 1)
    backend = StubBackend(log_text=_STEPS_EXCEEDED_LOG)
    backend.poll_return = "error"
    runner = JobRunner(conn, s, backend, io.StringIO())

    r = runner.resolve(
        _job(sig, job_kind="preopt", method="pm6", basis=""), "lbl", lambda: "com", wait=True
    )

    assert r.status == "complete"
    assert r.log_path is not None
    conn.close()


def test_resolve_preopt_real_error_stays_error(tmp_path):
    # SCF death / singular B-matrix (FormBX) Error-terminate WITHOUT the
    # steps-exceeded line → genuine failure, must stay 'error'.
    s = _settings(tmp_path)
    conn = db.open(s.db_path)
    sig = signature("c1ccccc1", "pm6", "", 0, 1)
    backend = StubBackend(
        log_text=" Convergence failure -- run terminated.\n Error termination ...\n"
    )
    backend.poll_return = "error"
    runner = JobRunner(conn, s, backend, io.StringIO())

    r = runner.resolve(
        _job(sig, job_kind="preopt", method="pm6", basis=""), "lbl", lambda: "com", wait=True
    )

    assert r.status == "error"
    conn.close()


# ── resolve_concurrent: lane-bounded parallel driver ─────────────────────────


class _LaneBackend:
    """Fake backend that records peak concurrency. Each job polls 'running' once,
    then takes `status_for(label)` (default 'complete'). submit/fetch bump an
    `active` counter so a test can assert the driver never exceeds `max_lanes`."""

    def __init__(self, *, status_for=None, log_provider=None):
        self.active = 0
        self.max_active = 0
        self.calls: list[tuple[str, str]] = []
        self._jobs: dict[str, str] = {}
        self._polls: dict[str, int] = {}
        self._label: dict[str, str] = {}
        self._status_for = status_for or (lambda label: "complete")
        self.log_provider = log_provider

    def submit(self, com_path, label):
        sid = f"stub-{label}"
        self.active += 1
        self.max_active = max(self.max_active, self.active)
        self._polls[sid] = 0
        self._label[sid] = label
        self.calls.append(("submit", label))
        return sid

    def poll(self, sid):
        self._polls[sid] += 1
        return "running" if self._polls[sid] <= 1 else self._status_for(self._label[sid])

    def fetch_log(self, sid, label, local_dir):
        self.active -= 1
        Path(local_dir).mkdir(parents=True, exist_ok=True)
        out = Path(local_dir) / f"{label}.log"
        out.write_text(self.log_provider(label) if self.log_provider else synth_log())
        return str(out)

    def cancel(self, sid):
        pass


def _spec(extras, label, **job_overrides):
    sig = signature("c1ccccc1", "b3lyp", "6-31g**", 0, 1, extras=extras)
    return JobSpec(job=_job(sig, **job_overrides), label=label, build_com=lambda: "com")


def test_resolve_concurrent_caps_lanes_and_returns_in_order(tmp_path):
    s = _settings(tmp_path)
    conn = db.open(s.db_path)
    energies = {"lbl0": -10.0, "lbl1": -20.0, "lbl2": -30.0}
    backend = _LaneBackend(log_provider=lambda label: synth_log(energy_ha=energies[label]))
    runner = JobRunner(conn, s, backend, io.StringIO())
    specs = [_spec(f"x{i}", f"lbl{i}") for i in range(3)]

    results = runner.resolve_concurrent(specs, max_lanes=2)

    assert backend.max_active == 2  # never more than 2 g16 in flight at once
    assert all(r.status == "complete" for r in results)
    # results stay in spec order — each carries its own log's energy.
    assert [round(r.result.energy_ev, 3) for r in results] == [
        round(e * HARTREE_TO_EV, 3) for e in (-10.0, -20.0, -30.0)
    ]
    conn.close()


def test_resolve_concurrent_serial_when_one_lane(tmp_path):
    s = _settings(tmp_path)
    conn = db.open(s.db_path)
    backend = _LaneBackend()
    runner = JobRunner(conn, s, backend, io.StringIO())
    specs = [_spec(f"x{i}", f"lbl{i}") for i in range(3)]

    results = runner.resolve_concurrent(specs, max_lanes=1)

    assert backend.max_active == 1  # max_lanes=1 → never overlaps (serial)
    assert all(r.status == "complete" for r in results)
    conn.close()


def test_resolve_concurrent_settles_cache_hits_without_submit(tmp_path):
    s = _settings(tmp_path)
    conn = db.open(s.db_path)
    sig_hit = signature("c1ccccc1", "b3lyp", "6-31g**", 0, 1, extras="hit")
    jid = db.insert_job(conn, _job(sig_hit, status="complete", log_path="/seeded.log"))
    db.insert_result(
        conn, db.Result(job_id=jid, homo_ev=-6.0, lumo_ev=-1.0, gap_ev=5.0, energy_ev=-4000.0)
    )
    backend = _LaneBackend()
    runner = JobRunner(conn, s, backend, io.StringIO())
    specs = [JobSpec(_job(sig_hit), "lbl_hit", lambda: "com"), _spec("miss", "lbl_miss")]

    results = runner.resolve_concurrent(specs, max_lanes=2)

    assert results[0].cached is True and results[0].result.energy_ev == -4000.0
    assert results[1].cached is False and results[1].status == "complete"
    # Only the cache miss was submitted.
    assert [c for c in backend.calls if c[0] == "submit"] == [("submit", "lbl_miss")]
    conn.close()


def test_resolve_concurrent_one_failure_doesnt_block_siblings(tmp_path):
    s = _settings(tmp_path)
    conn = db.open(s.db_path)
    backend = _LaneBackend(status_for=lambda label: "error" if label == "boom" else "complete")
    runner = JobRunner(conn, s, backend, io.StringIO())
    specs = [_spec("ok1", "ok1"), _spec("boom", "boom"), _spec("ok2", "ok2")]

    results = runner.resolve_concurrent(specs, max_lanes=2)

    assert [r.status for r in results] == ["complete", "error", "complete"]
    conn.close()


def test_resolve_concurrent_ctrl_c_detaches_whole_active_set(tmp_path):
    s = _settings(tmp_path)
    conn = db.open(s.db_path)

    class _Interrupt(_LaneBackend):
        def poll(self, sid):
            raise KeyboardInterrupt

    backend = _Interrupt()
    runner = JobRunner(conn, s, backend, io.StringIO())
    specs = [_spec(f"k{i}", f"k{i}") for i in range(2)]

    results = runner.resolve_concurrent(specs, max_lanes=2)

    assert all(r.status == "detached" for r in results)
    # Both rows stay 'running' with an ssh_jobid → reattachable next run (ADR 0005).
    running = conn.execute(
        "SELECT count(*) FROM jobs WHERE status='running' AND ssh_jobid IS NOT NULL"
    ).fetchone()[0]
    assert running == 2
    conn.close()
