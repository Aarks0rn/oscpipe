"""Tests for `oscpipe export` (flat one-row-per-candidate CSV)."""

from __future__ import annotations

import argparse
import csv
import io
import json

from oscpipe.cli.main import EXPORT_COLUMNS, run_export
from oscpipe.settings import Settings
from oscpipe.store import db


def _setup(tmp_path):
    s = Settings(
        backend="local",
        db_path=str(tmp_path / "results.db"),
        gaussian_nproc=2,
        gaussian_mem="2GB",
    )
    conn = db.open(s.db_path)
    return s, conn


def _job(signature, smiles, **overrides):
    base = dict(
        id=None,
        signature=signature,
        smiles=smiles,
        method="b3lyp",
        basis="6-31g*",
        charge=0,
        mult=1,
        job_kind="properties",
        status="complete",
        submitted_at="2026-06-11T00:00:00",
    )
    base.update(overrides)
    return db.Job(**base)


def _seed_candidate(conn, smiles, sig_prefix):
    jid = db.insert_job(conn, _job(f"{sig_prefix}-p", smiles))
    db.insert_result(conn, db.Result(job_id=jid, homo_ev=-5.2, lumo_ev=-2.9, gap_ev=2.3))
    tid = db.insert_job(conn, _job(f"{sig_prefix}-t", smiles, job_kind="tddft"))
    states = [{"n": 1, "energy_ev": 2.5, "wavelength_nm": 495.9, "f": 0.8}]
    db.insert_result(conn, db.Result(job_id=tid, spectra_json=json.dumps(states)))


def _read_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def test_export_writes_one_row_per_candidate(tmp_path):
    s, conn = _setup(tmp_path)
    _seed_candidate(conn, "c1ccccc1", "sig-a")
    _seed_candidate(conn, "c1ccsc1", "sig-b")

    out_csv = tmp_path / "out.csv"
    stdout = io.StringIO()
    rc = run_export(argparse.Namespace(csv=str(out_csv)), s, conn, stdout=stdout)
    assert rc == 0
    assert "export: 2 candidates" in stdout.getvalue()

    rows = _read_csv(out_csv)
    assert len(rows) == 2
    assert list(rows[0].keys()) == EXPORT_COLUMNS
    by_smiles = {r["smiles"]: r for r in rows}
    assert float(by_smiles["c1ccccc1"]["homo_ev"]) == -5.2
    assert float(by_smiles["c1ccccc1"]["lambda_max_nm"]) == 495.9


def test_export_blank_cells_for_missing_fields(tmp_path):
    s, conn = _setup(tmp_path)
    jid = db.insert_job(conn, _job("sig-only-p", "c1ccccc1"))
    db.insert_result(conn, db.Result(job_id=jid, homo_ev=-5.2, lumo_ev=-2.9, gap_ev=2.3))

    out_csv = tmp_path / "out.csv"
    rc = run_export(argparse.Namespace(csv=str(out_csv)), s, conn, stdout=io.StringIO())
    assert rc == 0
    rows = _read_csv(out_csv)
    assert len(rows) == 1
    assert rows[0]["lambda_max_nm"] == ""
    assert rows[0]["lambda_h_ev"] == ""


def test_export_empty_db_writes_header_only(tmp_path):
    s, conn = _setup(tmp_path)
    out_csv = tmp_path / "out.csv"
    rc = run_export(argparse.Namespace(csv=str(out_csv)), s, conn, stdout=io.StringIO())
    assert rc == 0
    with open(out_csv) as f:
        lines = f.read().strip().splitlines()
    assert len(lines) == 1
    assert lines[0].split(",") == EXPORT_COLUMNS
