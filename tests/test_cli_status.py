"""Tests for `oscpipe status` printing."""

from __future__ import annotations

import argparse

from oscpipe.cli.main import run_status
from oscpipe.settings import Settings
from oscpipe.store import db


def _settings(tmp_path):
    return Settings(db_path=str(tmp_path / "results.db"))


def _job(**over):
    base = dict(
        id=None,
        signature="sig",
        smiles="c1ccccc1",
        method="b3lyp",
        basis="6-31g*",
        charge=0,
        mult=1,
        job_kind="properties",
        status="pending",
        submitted_at="2026-05-21T00:00:00",
    )
    base.update(over)
    return db.Job(**base)


def test_status_empty_table(tmp_path, capsys):
    s = _settings(tmp_path)
    conn = db.open(s.db_path)
    rc = run_status(argparse.Namespace(), s, conn)
    assert rc == 0
    assert "(no jobs)" in capsys.readouterr().out


def test_status_lists_jobs(tmp_path, capsys):
    s = _settings(tmp_path)
    conn = db.open(s.db_path)
    j1 = db.insert_job(conn, _job(signature="s1", smiles="c1ccccc1", status="complete"))
    db.insert_result(
        conn,
        db.Result(job_id=j1, homo_ev=-6.7, lumo_ev=-1.0, gap_ev=5.7, dipole_debye=0.0),
    )
    db.insert_job(conn, _job(signature="s2", smiles="c1ccsc1", status="running"))

    run_status(argparse.Namespace(), s, conn)
    out = capsys.readouterr().out
    # Header + two rows.
    assert "status" in out
    assert "complete" in out
    assert "running" in out
    assert "c1ccccc1" in out
    assert "c1ccsc1" in out
    # HOMO column populated for the complete row, dashes for the running one.
    assert "-6.700" in out
    assert "—" in out
