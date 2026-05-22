"""Tests for `oscpipe reconcile`."""

from __future__ import annotations

import argparse

from _helpers import SYNTH_LOG, StubBackend, seed_running_job

from oscpipe.cli.main import run_reconcile
from oscpipe.settings import Settings
from oscpipe.store import db


def _settings(tmp_path):
    return Settings(
        backend="local",
        db_path=str(tmp_path / "results.db"),
        remote_work_dir="/work",
    )


def test_reconcile_empty(tmp_path, capsys):
    s = _settings(tmp_path)
    conn = db.open(s.db_path)
    rc = run_reconcile(argparse.Namespace(), s, StubBackend(), conn)
    assert rc == 0
    assert "nothing to do" in capsys.readouterr().out


def test_reconcile_promotes_completed_job(tmp_path):
    s = _settings(tmp_path)
    conn = db.open(s.db_path)
    jid = seed_running_job(conn)
    backend = StubBackend(log_text=SYNTH_LOG)
    backend.poll_return = "complete"

    rc = run_reconcile(argparse.Namespace(), s, backend, conn)
    assert rc == 0
    row = conn.execute("SELECT status FROM jobs WHERE id = ?", (jid,)).fetchone()
    assert row["status"] == "complete"


def test_reconcile_marks_unknown_jobs_lost(tmp_path, capsys):
    s = _settings(tmp_path)
    conn = db.open(s.db_path)
    jid = seed_running_job(conn)
    backend = StubBackend()
    backend.poll_return = "unknown"

    run_reconcile(argparse.Namespace(), s, backend, conn)
    row = conn.execute("SELECT * FROM jobs WHERE id = ?", (jid,)).fetchone()
    assert row["status"] == "error"
    assert "lost" in (row["error_msg"] or "")
    assert "1 marked lost" in capsys.readouterr().out


def test_reconcile_handles_job_without_ssh_jobid(tmp_path):
    s = _settings(tmp_path)
    conn = db.open(s.db_path)
    jid = db.insert_job(
        conn,
        db.Job(
            id=None,
            signature="x",
            smiles="c1ccccc1",
            method="b3lyp",
            basis="6-31g*",
            charge=0,
            mult=1,
            job_kind="properties",
            status="pending",
            submitted_at="2026-05-21T00:00:00",
        ),
    )
    run_reconcile(argparse.Namespace(), s, StubBackend(), conn)
    row = conn.execute("SELECT * FROM jobs WHERE id = ?", (jid,)).fetchone()
    assert row["status"] == "error"
    assert "no ssh_jobid" in (row["error_msg"] or "")
