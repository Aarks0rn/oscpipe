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


def test_status_shows_charge_column(tmp_path, capsys):
    """status shows each job's charge so a charged sub-job's orbitals (e.g. a λ_h
    cation_opt) are not misread as neutral properties."""
    s = _settings(tmp_path)
    conn = db.open(s.db_path)
    jid = db.insert_job(conn, _job(signature="cat", charge=1, mult=2, status="complete"))
    db.insert_result(conn, db.Result(job_id=jid, homo_ev=-9.9, lumo_ev=-5.0, gap_ev=4.9))

    run_status(argparse.Namespace(), s, conn)
    out = capsys.readouterr().out
    assert "chg" in out  # charge column header present
    data = next(ln for ln in out.splitlines() if "c1ccccc1" in ln)
    fields = data.split()  # id  status  method/basis  chg  HOMO  gap  smiles
    assert fields[3] == "1"
    assert fields[4] == "-9.900"


class _StatusBackend:
    """Minimal backend for the liveness column: poll returns a fixed PID state."""

    def __init__(self, status="running"):
        self._status = status
        self._jobs: dict[str, str] = {}

    def poll(self, sid):
        return self._status


def test_status_running_filter_hides_terminal(tmp_path, capsys):
    s = _settings(tmp_path)
    conn = db.open(s.db_path)
    db.insert_job(conn, _job(signature="done", smiles="c1ccccc1", status="complete"))
    db.insert_job(conn, _job(signature="run", smiles="c1ccsc1", status="running", ssh_jobid="111"))
    db.insert_job(conn, _job(signature="pend", smiles="c1ccoc1", status="pending"))

    run_status(argparse.Namespace(running=True), s, conn)
    out = capsys.readouterr().out
    assert "c1ccsc1" in out and "c1ccoc1" in out  # running + pending shown
    assert "c1ccccc1" not in out  # complete hidden


def test_status_liveness_column_only_with_backend(tmp_path, capsys):
    s = _settings(tmp_path)
    conn = db.open(s.db_path)
    db.insert_job(
        conn,
        _job(
            signature="run",
            smiles="c1ccsc1",
            status="running",
            started_at="2026-05-21T00:00:00",
            ssh_jobid="111",
        ),
    )

    run_status(argparse.Namespace(), s, conn)  # no backend → no live column
    assert "live" not in capsys.readouterr().out

    run_status(argparse.Namespace(), s, conn, backend=_StatusBackend(status="running"))
    out = capsys.readouterr().out
    assert "live" in out  # live column header present
    assert "ALIVE" in out  # polled PID state rendered
