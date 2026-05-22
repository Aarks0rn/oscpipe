"""Tests for `oscpipe fetch`.

A StubBackend is injected so the lifecycle can be exercised deterministically
without spawning subprocesses or talking to a workstation. parse_properties
runs against a tiny synthetic Gaussian log written into local_dir on the fly.
"""

from __future__ import annotations

import argparse

from _helpers import EMPTY_LOG, SYNTH_LOG, StubBackend, seed_running_job

from oscpipe.cli.main import run_fetch
from oscpipe.settings import Settings
from oscpipe.store import db


def _settings(tmp_path):
    return Settings(
        backend="local",
        db_path=str(tmp_path / "results.db"),
        gaussian_nproc=2,
        gaussian_mem="2GB",
    )


def _ns(job_id=None):
    return argparse.Namespace(job_id=job_id)


# ── empty case ─────────────────────────────────────────────────────────────


def test_fetch_no_running_jobs(tmp_path, capsys):
    s = _settings(tmp_path)
    conn = db.open(s.db_path)
    rc = run_fetch(_ns(), s, StubBackend(), conn)
    assert rc == 0
    assert "no running jobs" in capsys.readouterr().out


# ── happy path ─────────────────────────────────────────────────────────────


def test_fetch_completes_and_persists_result(tmp_path, capsys):
    s = _settings(tmp_path)
    conn = db.open(s.db_path)
    jid = seed_running_job(conn)
    backend = StubBackend(log_text=SYNTH_LOG)
    backend.poll_return = "complete"

    rc = run_fetch(_ns(), s, backend, conn)
    assert rc == 0

    row = conn.execute("SELECT * FROM jobs WHERE id = ?", (jid,)).fetchone()
    assert row["status"] == "complete"
    assert row["completed_at"] is not None
    assert row["log_path"].endswith(".log")

    result = conn.execute("SELECT * FROM results WHERE job_id = ?", (jid,)).fetchone()
    assert result is not None
    # HOMO = -0.5 Ha = -13.6057 eV, LUMO = 0.1 Ha = 2.72114 eV
    assert abs(result["homo_ev"] - (-13.6057)) < 1e-2
    assert abs(result["lumo_ev"] - 2.7211) < 1e-2
    assert abs(result["gap_ev"] - 16.3268) < 1e-2

    out = capsys.readouterr().out
    assert "complete: job=" in out
    assert "1 complete, 0 error, 0 still running" in out


# ── parse failure ──────────────────────────────────────────────────────────


def test_fetch_parse_failure_marks_error(tmp_path, capsys):
    s = _settings(tmp_path)
    conn = db.open(s.db_path)
    jid = seed_running_job(conn)
    backend = StubBackend(log_text=EMPTY_LOG)
    backend.poll_return = "complete"

    rc = run_fetch(_ns(), s, backend, conn)
    assert rc == 0

    row = conn.execute("SELECT * FROM jobs WHERE id = ?", (jid,)).fetchone()
    assert row["status"] == "error"
    assert "orbital" in (row["error_msg"] or "")


# ── still running ──────────────────────────────────────────────────────────


def test_fetch_still_running_keeps_status(tmp_path, capsys):
    s = _settings(tmp_path)
    conn = db.open(s.db_path)
    jid = seed_running_job(conn)
    backend = StubBackend()
    backend.poll_return = "running"
    run_fetch(_ns(), s, backend, conn)
    row = conn.execute("SELECT * FROM jobs WHERE id = ?", (jid,)).fetchone()
    assert row["status"] == "running"
    out = capsys.readouterr().out
    assert "1 still running" in out


def test_fetch_backend_error_marks_error(tmp_path):
    s = _settings(tmp_path)
    conn = db.open(s.db_path)
    jid = seed_running_job(conn)
    backend = StubBackend()
    backend.poll_return = "error"
    run_fetch(_ns(), s, backend, conn)
    row = conn.execute("SELECT * FROM jobs WHERE id = ?", (jid,)).fetchone()
    assert row["status"] == "error"
    assert "backend reported error" in (row["error_msg"] or "")


# ── --job-id filter ────────────────────────────────────────────────────────


def test_fetch_filters_by_job_id(tmp_path):
    s = _settings(tmp_path)
    conn = db.open(s.db_path)
    j1 = seed_running_job(conn, sig="sig-a", ssh_jobid="j-a")
    j2 = seed_running_job(conn, sig="sig-b", ssh_jobid="j-b")
    backend = StubBackend(log_text=SYNTH_LOG)
    backend.poll_return = "complete"

    run_fetch(_ns(job_id=j1), s, backend, conn)

    r1 = conn.execute("SELECT status FROM jobs WHERE id = ?", (j1,)).fetchone()
    r2 = conn.execute("SELECT status FROM jobs WHERE id = ?", (j2,)).fetchone()
    assert r1["status"] == "complete"
    assert r2["status"] == "running"
    # Backend.poll was called exactly once (for j1 only).
    assert sum(1 for c in backend.calls if c[0] == "poll") == 1


# ── ssh rehydration code path ──────────────────────────────────────────────


def test_fetch_rehydrates_ssh_backend_jobs_map(tmp_path):
    s = Settings(
        backend="local",
        db_path=str(tmp_path / "results.db"),
        remote_work_dir="/home/user/work",
    )
    conn = db.open(s.db_path)
    seed_running_job(conn, ssh_jobid="555")
    backend = StubBackend()
    backend.poll_return = "running"  # no work happens, just rehydrate
    run_fetch(_ns(), s, backend, conn)
    assert "555" in backend._jobs
    assert backend._jobs["555"].endswith(".log")
    assert backend._jobs["555"].startswith("/home/user/work/")
