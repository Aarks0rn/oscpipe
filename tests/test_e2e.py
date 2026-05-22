"""End-to-end smoke: SMILES → backend → DB row.

All tests use StubBackend so the suite runs offline.
Real workstation tests (SshBackend + qsub) live in test_real_workstation.py.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from _helpers import StubBackend

from oscpipe.cli.main import run_submit
from oscpipe.settings import Settings
from oscpipe.store import db

_FIXTURES = Path(__file__).parent / "fixtures"

_MINIMAL_LOG = " Normal termination of Gaussian 16 at ...\n"


def _args(smiles: str):
    return argparse.Namespace(smiles=smiles, method="b3lyp", basis="6-31g*", charge=0, mult=1)


def _setup(tmp_path, log_text=_MINIMAL_LOG, poll_return="complete"):
    s = Settings(
        backend="local",
        db_path=str(tmp_path / "results.db"),
        gaussian_nproc=2,
        gaussian_mem="2GB",
    )
    backend = StubBackend(log_text=log_text)
    backend.poll_return = poll_return
    conn = db.open(s.db_path)
    return s, backend, conn


def test_e2e_submit_lifecycle_records_job(tmp_path):
    """Minimal log (no orbital data) → parse fails → status=error, rc != 0."""
    s, backend, conn = _setup(tmp_path, log_text=_MINIMAL_LOG)
    rc = run_submit(_args("c1ccccc1"), s, backend, conn)
    assert rc != 0
    row = conn.execute("SELECT * FROM jobs").fetchone()
    assert row["smiles"] == "c1ccccc1"
    assert row["status"] == "error"
    assert row["error_msg"]


def test_e2e_cache_hit_short_circuits(tmp_path, capsys):
    s, backend, conn = _setup(tmp_path)
    from oscpipe.chem.smiles import canonicalise
    from oscpipe.store.cache import signature

    canon, _ = canonicalise("c1ccccc1")
    sig = signature(canon, "b3lyp", "6-31g*", 0, 1)
    jid = db.insert_job(
        conn,
        db.Job(
            id=None,
            signature=sig,
            smiles=canon,
            method="b3lyp",
            basis="6-31g*",
            charge=0,
            mult=1,
            job_kind="properties",
            status="complete",
            submitted_at="2026-05-22T00:00:00",
        ),
    )
    db.insert_result(conn, db.Result(job_id=jid, homo_ev=-6.7, lumo_ev=-1.0, gap_ev=5.7))

    rc = run_submit(_args("c1ccccc1"), s, backend, conn)
    assert rc == 0
    assert conn.execute("SELECT count(*) FROM jobs").fetchone()[0] == 1
    assert "cache hit" in capsys.readouterr().out


def test_e2e_warnings_are_printed(tmp_path, capsys):
    """Fluorene → canonicalise emits an sp3-in-aromatic warning."""
    s, backend, conn = _setup(tmp_path)
    run_submit(_args("c1ccc2c(c1)Cc1ccccc1-2"), s, backend, conn)
    out = capsys.readouterr().out
    assert "warning:" in out
    assert "sp3 carbon" in out


def test_e2e_happy_path_real_fixture(tmp_path, capsys):
    """Full SMILES → parse → DB using real h2.log (HF/STO-3G from workstation)."""
    h2_log = _FIXTURES / "h2.log"
    if not h2_log.exists():
        import pytest

        pytest.skip("tests/fixtures/h2.log not available")

    s, backend, conn = _setup(tmp_path, log_text=h2_log.read_text())
    rc = run_submit(_args("[H][H]"), s, backend, conn)
    assert rc == 0

    row = conn.execute("SELECT * FROM jobs ORDER BY id").fetchone()
    assert row["status"] == "complete"

    result = conn.execute(
        "SELECT homo_ev, lumo_ev, gap_ev, dipole_debye FROM results WHERE job_id = ?",
        (row["id"],),
    ).fetchone()
    assert abs(result["homo_ev"] - (-15.732)) < 1e-2
    assert abs(result["lumo_ev"] - 18.235) < 1e-2
    assert result["gap_ev"] > 0
    assert abs(result["dipole_debye"] - 0.0) < 1e-2

    out = capsys.readouterr().out
    assert "complete: job=1" in out
    assert "HOMO=" in out
