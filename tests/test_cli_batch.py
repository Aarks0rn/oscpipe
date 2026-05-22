"""Tests for `oscpipe batch`."""

from __future__ import annotations

import argparse
import io

from _helpers import StubBackend

from oscpipe.cli.main import run_batch
from oscpipe.settings import Settings
from oscpipe.store import db


def _args(file: str):
    return argparse.Namespace(
        file=file,
        method="b3lyp",
        basis="6-31g*",
        charge=0,
        mult=1,
    )


def _setup(tmp_path):
    s = Settings(
        backend="local",
        db_path=str(tmp_path / "results.db"),
        gaussian_nproc=2,
        gaussian_mem="2GB",
    )
    backend = StubBackend()
    backend.poll_return = "complete"
    conn = db.open(s.db_path)
    return s, backend, conn


def test_batch_reads_csv_with_header(tmp_path):
    s, backend, conn = _setup(tmp_path)
    csv_file = tmp_path / "list.csv"
    csv_file.write_text("smiles\nc1ccccc1\nc1ccsc1\n")
    rc = run_batch(_args(str(csv_file)), s, backend, conn)
    assert rc == 0
    n = conn.execute("SELECT count(*) FROM jobs").fetchone()[0]
    assert n == 2


def test_batch_groups_jobs_under_one_workflow(tmp_path):
    s, backend, conn = _setup(tmp_path)
    csv_file = tmp_path / "list.csv"
    csv_file.write_text("smiles\nc1ccccc1\nc1ccsc1\n")
    run_batch(_args(str(csv_file)), s, backend, conn)

    wf = conn.execute("SELECT * FROM workflows WHERE kind='batch'").fetchone()
    assert wf is not None
    assert wf["status"] == "complete"

    wf_ids = {
        r[0] for r in conn.execute("SELECT workflow_id FROM jobs WHERE workflow_id IS NOT NULL")
    }
    assert wf_ids == {wf["id"]}


def test_batch_workflow_summary_json(tmp_path):
    import json

    s, backend, conn = _setup(tmp_path)
    csv_file = tmp_path / "list.csv"
    csv_file.write_text("smiles\nc1ccccc1\nc1ccsc1\nCC\n")
    run_batch(_args(str(csv_file)), s, backend, conn)

    wf = conn.execute("SELECT summary_json FROM workflows WHERE kind='batch'").fetchone()
    summary = json.loads(wf["summary_json"])
    assert summary["total"] == 3
    assert summary["ok"] == 3
    assert summary["failed"] == 0


def test_batch_reads_stdin_plain_list(tmp_path):
    s, backend, conn = _setup(tmp_path)
    stdin = io.StringIO("c1ccccc1\nc1ccsc1\n\n")
    run_batch(_args("-"), s, backend, conn, stdin=stdin)
    n = conn.execute("SELECT count(*) FROM jobs").fetchone()[0]
    assert n == 2


def test_batch_reads_plain_file_without_header(tmp_path):
    s, backend, conn = _setup(tmp_path)
    f = tmp_path / "list.txt"
    f.write_text("c1ccccc1\nc1ccsc1\n")
    run_batch(_args(str(f)), s, backend, conn)
    n = conn.execute("SELECT count(*) FROM jobs").fetchone()[0]
    assert n == 2


def test_batch_empty_returns_error(tmp_path, capsys):
    s, backend, conn = _setup(tmp_path)
    f = tmp_path / "empty.txt"
    f.write_text("")
    rc = run_batch(_args(str(f)), s, backend, conn)
    assert rc == 1
    assert "no SMILES" in capsys.readouterr().out
