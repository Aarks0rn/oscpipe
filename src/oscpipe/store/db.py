"""Raw sqlite3 store. Two tables (jobs, results) + workflows + view.
Schema lives in schema.sql — open() applies it idempotently.
"""

from __future__ import annotations

import json
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
    # Concurrency posture. WAL lets the Streamlit dashboard (and any second oscpipe
    # process) read while the CLI writes without "database is locked"; busy_timeout
    # makes a writer wait for a competing write rather than erroring immediately.
    # WAL is persistent on the file; busy_timeout is per-connection. Needs a local
    # filesystem (results.db is local) — do not move the DB to a network mount.
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA busy_timeout=5000")
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


def find_inflight_by_signature(conn: sqlite3.Connection, signature: str):
    """Return the most recent pending/running job row with this signature that
    has a recorded ssh_jobid, or None.

    Used to reattach to an already-queued job instead of submitting a duplicate
    — the main source of orphan UGE jobs is a detached run (Ctrl-C / disconnect)
    that leaves a 'running' row + live remote job, then a re-run resubmitting
    because the cache only matches 'complete'.
    """
    return conn.execute(
        "SELECT * FROM jobs "
        "WHERE signature = ? AND status IN ('pending', 'running') "
        "AND ssh_jobid IS NOT NULL ORDER BY id DESC LIMIT 1",
        (signature,),
    ).fetchone()


# ── read API for external consumers (project scripts, `oscpipe export`) ─────
# The one home for "what shape do results take when they leave the store" —
# scripts use these instead of writing raw SQL against the schema.


def list_complete_oligomer_sweeps(conn: sqlite3.Connection) -> list[dict]:
    """Parsed summary of every complete oligomer_sweep workflow, oldest first.

    Each dict is the workflow's summary_json ("repeat_unit", "max_n", "per_n",
    "extrapolated", ...) plus "workflow_id" and "smiles" from the row.  Rows
    with no summary are skipped.  Oldest-first iteration gives callers keying
    by repeat unit newest-wins semantics.
    """
    out: list[dict] = []
    for row in conn.execute(
        "SELECT id, smiles, summary_json FROM workflows "
        "WHERE kind = 'oligomer_sweep' AND status = 'complete' "
        "AND summary_json IS NOT NULL ORDER BY id"
    ):
        summary = json.loads(row["summary_json"])
        summary["workflow_id"] = row["id"]
        summary["smiles"] = row["smiles"]
        out.append(summary)
    return out


def candidate_summary(conn: sqlite3.Connection, smiles: str) -> dict:
    """Latest complete results for one canonical SMILES, joined across the
    three sources: properties job (HOMO/LUMO/gap), tddft job (brightest state),
    lambda_h workflow (λ_h / J / Marcus rate).  Missing pieces are None.

    No physics here — this only reshapes stored rows.
    """
    out: dict = {
        "smiles": smiles,
        "homo_ev": None,
        "lumo_ev": None,
        "gap_ev": None,
        "lambda_max_nm": None,
        "f_osc": None,
        "lambda_h_ev": None,
        "j_hole_ev": None,
        "marcus_rate_s1": None,
    }
    # charge=0 keeps cation_opt rows (also job_kind='properties') out.
    row = conn.execute(
        "SELECT r.homo_ev, r.lumo_ev, r.gap_ev FROM jobs j "
        "JOIN results r ON r.job_id = j.id "
        "WHERE j.smiles = ? AND j.job_kind = 'properties' AND j.charge = 0 "
        "AND j.status = 'complete' ORDER BY j.id DESC LIMIT 1",
        (smiles,),
    ).fetchone()
    if row is not None:
        out["homo_ev"] = row["homo_ev"]
        out["lumo_ev"] = row["lumo_ev"]
        out["gap_ev"] = row["gap_ev"]

    row = conn.execute(
        "SELECT r.spectra_json FROM jobs j JOIN results r ON r.job_id = j.id "
        "WHERE j.smiles = ? AND j.job_kind = 'tddft' AND j.status = 'complete' "
        "AND r.spectra_json IS NOT NULL ORDER BY j.id DESC LIMIT 1",
        (smiles,),
    ).fetchone()
    if row is not None:
        states = json.loads(row["spectra_json"])
        if states:
            bright = max(states, key=lambda s: s["f"])
            out["lambda_max_nm"] = bright["wavelength_nm"]
            out["f_osc"] = bright["f"]

    row = conn.execute(
        "SELECT summary_json FROM workflows "
        "WHERE kind = 'lambda_h' AND smiles = ? AND status = 'complete' "
        "AND summary_json IS NOT NULL ORDER BY id DESC LIMIT 1",
        (smiles,),
    ).fetchone()
    if row is not None:
        summary = json.loads(row["summary_json"])
        out["lambda_h_ev"] = summary.get("lambda_h_ev")
        out["j_hole_ev"] = summary.get("j_hole_ev")
        out["marcus_rate_s1"] = summary.get("marcus_rate_s1")
    return out


def list_candidate_smiles(conn: sqlite3.Connection) -> list[str]:
    """Distinct SMILES with any complete charge-0 properties job, complete
    tddft job, or complete lambda_h workflow, ordered by SMILES.

    Oligomer n-mer SMILES appear too — they are real molecules with real
    results; consumers filter if they only want repeat units.
    """
    rows = conn.execute(
        "SELECT smiles FROM jobs "
        "WHERE job_kind = 'properties' AND charge = 0 AND status = 'complete' "
        "UNION "
        "SELECT smiles FROM jobs WHERE job_kind = 'tddft' AND status = 'complete' "
        "UNION "
        "SELECT smiles FROM workflows WHERE kind = 'lambda_h' AND status = 'complete' "
        "ORDER BY smiles"
    )
    return [r["smiles"] for r in rows]
