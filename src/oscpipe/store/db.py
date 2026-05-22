"""Raw sqlite3 store. Two tables (jobs, results) + workflows + view.
Schema lives in schema.sql — open() applies it idempotently.
"""

from __future__ import annotations

import sqlite3
from dataclasses import dataclass
from pathlib import Path

SCHEMA = Path(__file__).with_name("schema.sql")


@dataclass
class Job:
    id: int | None
    signature: str
    smiles: str
    method: str
    basis: str
    charge: int
    mult: int
    job_kind: str
    status: str
    submitted_at: str
    started_at: str | None = None
    completed_at: str | None = None
    log_path: str | None = None
    remote_log_path: str | None = None
    ssh_jobid: str | None = None
    error_msg: str | None = None
    workflow_id: int | None = None
    notes: str | None = None


@dataclass
class Result:
    job_id: int
    homo_ev: float | None = None
    lumo_ev: float | None = None
    gap_ev: float | None = None
    dipole_debye: float | None = None
    energy_ev: float | None = None
    atoms_xyz: str | None = None
    spectra_json: str | None = None


def open(db_path: str) -> sqlite3.Connection:
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    conn.executescript(SCHEMA.read_text())
    return conn


_JOB_COLUMNS = (
    "signature",
    "smiles",
    "method",
    "basis",
    "charge",
    "mult",
    "job_kind",
    "status",
    "submitted_at",
    "started_at",
    "completed_at",
    "log_path",
    "remote_log_path",
    "ssh_jobid",
    "error_msg",
    "workflow_id",
    "notes",
)

_RESULT_COLUMNS = (
    "job_id",
    "homo_ev",
    "lumo_ev",
    "gap_ev",
    "dipole_debye",
    "energy_ev",
    "atoms_xyz",
    "spectra_json",
)


def insert_job(conn: sqlite3.Connection, job: Job) -> int:
    cols = ", ".join(_JOB_COLUMNS)
    placeholders = ", ".join("?" for _ in _JOB_COLUMNS)
    values = tuple(getattr(job, c) for c in _JOB_COLUMNS)
    with conn:
        cur = conn.execute(f"INSERT INTO jobs ({cols}) VALUES ({placeholders})", values)
    return cur.lastrowid


def update_job_status(conn: sqlite3.Connection, job_id: int, status: str, **fields) -> None:
    unknown = set(fields) - set(_JOB_COLUMNS)
    if unknown:
        raise ValueError(f"unknown job columns: {sorted(unknown)}")
    assignments = ["status = ?"] + [f"{k} = ?" for k in fields]
    values = [status, *fields.values(), job_id]
    with conn:
        conn.execute(f"UPDATE jobs SET {', '.join(assignments)} WHERE id = ?", values)


def insert_result(conn: sqlite3.Connection, result: Result) -> None:
    cols = ", ".join(_RESULT_COLUMNS)
    placeholders = ", ".join("?" for _ in _RESULT_COLUMNS)
    values = tuple(getattr(result, c) for c in _RESULT_COLUMNS)
    with conn:
        conn.execute(f"INSERT INTO results ({cols}) VALUES ({placeholders})", values)


_WORKFLOW_COLUMNS = ("kind", "smiles", "created_at", "status", "summary_json")


def insert_workflow(
    conn: sqlite3.Connection,
    kind: str,
    smiles: str,
    created_at: str,
    status: str = "pending",
) -> int:
    with conn:
        cur = conn.execute(
            "INSERT INTO workflows (kind, smiles, created_at, status) VALUES (?, ?, ?, ?)",
            (kind, smiles, created_at, status),
        )
    return cur.lastrowid


def update_workflow(conn: sqlite3.Connection, workflow_id: int, status: str, **fields) -> None:
    unknown = set(fields) - set(_WORKFLOW_COLUMNS)
    if unknown:
        raise ValueError(f"unknown workflow columns: {sorted(unknown)}")
    assignments = ["status = ?"] + [f"{k} = ?" for k in fields]
    values = [status, *fields.values(), workflow_id]
    with conn:
        conn.execute(f"UPDATE workflows SET {', '.join(assignments)} WHERE id = ?", values)


def find_complete_by_signature(conn: sqlite3.Connection, signature: str):
    """Return the joined jobs+results row for a complete job, or None."""
    return conn.execute(
        "SELECT * FROM v_jobs_with_results "
        "WHERE signature = ? AND status = 'complete' "
        "ORDER BY id DESC LIMIT 1",
        (signature,),
    ).fetchone()
