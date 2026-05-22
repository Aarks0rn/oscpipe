"""Store tests — real SQLite in tmp_path, no mocks."""

from oscpipe.store import db


def _job(**overrides):
    base = dict(
        id=None,
        signature="sig-aaa",
        smiles="c1ccccc1",
        method="b3lyp",
        basis="6-31g*",
        charge=0,
        mult=1,
        job_kind="properties",
        status="pending",
        submitted_at="2026-05-21T10:00:00",
    )
    base.update(overrides)
    return db.Job(**base)


# ── schema + signature (pre-existing) ──────────────────────────────────────


def test_open_creates_schema(tmp_path):
    conn = db.open(str(tmp_path / "results.db"))
    tables = {r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    assert {"jobs", "results", "workflows"}.issubset(tables)
    conn.close()


def test_signature_is_stable():
    from oscpipe.store.cache import signature

    a = signature("c1ccsc1", "b3lyp", "6-31g*", 0, 1)
    b = signature("c1ccsc1", "B3LYP", "6-31G*", 0, 1)
    assert a == b  # case-insensitive on method/basis


# ── insert_job ─────────────────────────────────────────────────────────────


def test_insert_job_returns_rowid_and_persists(tmp_path):
    conn = db.open(str(tmp_path / "r.db"))
    jid = db.insert_job(conn, _job())
    assert isinstance(jid, int) and jid > 0
    row = conn.execute("SELECT * FROM jobs WHERE id = ?", (jid,)).fetchone()
    assert row["signature"] == "sig-aaa"
    assert row["status"] == "pending"
    conn.close()


# ── update_job_status ──────────────────────────────────────────────────────


def test_update_job_status_changes_status(tmp_path):
    conn = db.open(str(tmp_path / "r.db"))
    jid = db.insert_job(conn, _job())
    db.update_job_status(conn, jid, "running", started_at="2026-05-21T10:01:00")
    row = conn.execute("SELECT * FROM jobs WHERE id = ?", (jid,)).fetchone()
    assert row["status"] == "running"
    assert row["started_at"] == "2026-05-21T10:01:00"
    conn.close()


def test_update_job_status_ignores_unknown_field(tmp_path):
    conn = db.open(str(tmp_path / "r.db"))
    jid = db.insert_job(conn, _job())
    # Unknown columns should be rejected, not silently swallowed.
    import pytest

    with pytest.raises(ValueError):
        db.update_job_status(conn, jid, "running", not_a_real_column="x")
    conn.close()


# ── insert_result ──────────────────────────────────────────────────────────


def test_insert_result_links_to_job(tmp_path):
    conn = db.open(str(tmp_path / "r.db"))
    jid = db.insert_job(conn, _job())
    db.insert_result(
        conn,
        db.Result(
            job_id=jid,
            homo_ev=-5.4,
            lumo_ev=-2.1,
            gap_ev=3.3,
            dipole_debye=0.0,
            energy_ev=-1000.0,
        ),
    )
    row = conn.execute("SELECT * FROM results WHERE job_id = ?", (jid,)).fetchone()
    assert row["homo_ev"] == -5.4
    assert row["gap_ev"] == 3.3
    conn.close()


# ── find_complete_by_signature ─────────────────────────────────────────────


def test_find_complete_by_signature_returns_none_when_missing(tmp_path):
    conn = db.open(str(tmp_path / "r.db"))
    assert db.find_complete_by_signature(conn, "nope") is None
    conn.close()


def test_find_complete_by_signature_ignores_incomplete(tmp_path):
    conn = db.open(str(tmp_path / "r.db"))
    db.insert_job(conn, _job(signature="sig-x", status="pending"))
    assert db.find_complete_by_signature(conn, "sig-x") is None
    conn.close()


def test_find_complete_by_signature_returns_joined_row(tmp_path):
    conn = db.open(str(tmp_path / "r.db"))
    jid = db.insert_job(conn, _job(signature="sig-y", status="complete"))
    db.insert_result(conn, db.Result(job_id=jid, homo_ev=-5.0, lumo_ev=-2.0, gap_ev=3.0))
    row = db.find_complete_by_signature(conn, "sig-y")
    assert row is not None
    assert row["homo_ev"] == -5.0
    assert row["gap_ev"] == 3.0
    assert row["status"] == "complete"
    conn.close()
