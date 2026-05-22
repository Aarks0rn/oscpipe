"""Tests for `oscpipe preflight` workstation health-check subcommand."""

from __future__ import annotations

import argparse

from _helpers import StubBackend

from oscpipe.cli.main import run_preflight
from oscpipe.settings import Settings
from oscpipe.store import db


def _setup(tmp_path):
    s = Settings(backend="local", db_path=str(tmp_path / "results.db"))
    conn = db.open(s.db_path)
    return s, conn


class _PreflightBackend(StubBackend):
    def __init__(self, checks: list[tuple[str, bool, str]]):
        super().__init__()
        self._checks = checks

    def preflight(self) -> list[tuple[str, bool, str]]:
        return self._checks


_ALL_PASS = [
    ("g16", True, "/opt/g16/g16"),
    ("qstat", True, "Following jobs do not exist: 0"),
    ("scratch", True, "/scratch/user"),
]


def test_preflight_all_pass(tmp_path, capsys):
    s, conn = _setup(tmp_path)
    rc = run_preflight(argparse.Namespace(), s, _PreflightBackend(_ALL_PASS), conn)
    assert rc == 0
    out = capsys.readouterr().out
    assert "db: ok" in out
    assert out.count("[ok]") == 3
    assert "FAIL" not in out


def test_preflight_g16_missing(tmp_path, capsys):
    s, conn = _setup(tmp_path)
    checks = [
        ("g16", False, "not found in PATH"),
        ("qstat", True, "Following jobs do not exist: 0"),
        ("scratch", True, "/scratch/user"),
    ]
    rc = run_preflight(argparse.Namespace(), s, _PreflightBackend(checks), conn)
    assert rc == 1
    out = capsys.readouterr().out
    assert "g16: [FAIL]" in out
    assert "not found in PATH" in out


def test_preflight_scratch_not_writable(tmp_path, capsys):
    s, conn = _setup(tmp_path)
    checks = [
        ("g16", True, "/opt/g16/g16"),
        ("qstat", True, "Following jobs do not exist: 0"),
        ("scratch", False, "not writable: Permission denied"),
    ]
    rc = run_preflight(argparse.Namespace(), s, _PreflightBackend(checks), conn)
    assert rc == 1
    out = capsys.readouterr().out
    assert "scratch: [FAIL]" in out


def test_preflight_multiple_failures(tmp_path, capsys):
    s, conn = _setup(tmp_path)
    checks = [
        ("g16", False, "not found in PATH"),
        ("qstat", False, "qstat not found"),
        ("scratch", False, "not writable: Permission denied"),
    ]
    rc = run_preflight(argparse.Namespace(), s, _PreflightBackend(checks), conn)
    assert rc == 1
    out = capsys.readouterr().out
    assert out.count("[FAIL]") == 3


def test_preflight_no_preflight_method(tmp_path, capsys):
    """Backend without preflight() (plain StubBackend) → graceful skip, rc=0."""
    s, conn = _setup(tmp_path)
    rc = run_preflight(argparse.Namespace(), s, StubBackend(), conn)
    assert rc == 0
    out = capsys.readouterr().out
    assert "db: ok" in out
    assert "no preflight" in out


def test_preflight_db_reports_job_count(tmp_path, capsys):
    s, conn = _setup(tmp_path)
    # Seed two jobs so count > 0.
    from oscpipe.store import db as _db
    from oscpipe.chem.smiles import canonicalise

    for smi in ("c1ccccc1", "CC"):
        canon, _ = canonicalise(smi)
        _db.insert_job(
            conn,
            _db.Job(
                id=None,
                signature="sig",
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

    rc = run_preflight(argparse.Namespace(), s, StubBackend(), conn)
    assert rc == 0
    assert "db: ok (2 jobs)" in capsys.readouterr().out
