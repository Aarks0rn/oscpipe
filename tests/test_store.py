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


def test_open_enables_wal_and_busy_timeout(tmp_path):
    conn = db.open(str(tmp_path / "results.db"))
    assert conn.execute("PRAGMA journal_mode").fetchone()[0] == "wal"
    assert conn.execute("PRAGMA busy_timeout").fetchone()[0] == 5000
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


# ── read API: list_complete_oligomer_sweeps / candidate_summary / list_candidate_smiles


def _seed_workflow(conn, kind, smiles, status="complete", summary=None):
    import json

    wf_id = db.insert_workflow(conn, kind, smiles, "2026-06-11T00:00:00", status=status)
    if summary is not None:
        db.update_workflow(conn, wf_id, status, summary_json=json.dumps(summary))
    return wf_id


def test_list_complete_oligomer_sweeps_filters_and_parses(tmp_path):
    conn = db.open(str(tmp_path / "r.db"))
    ok_id = _seed_workflow(
        conn, "oligomer_sweep", "ru-1", summary={"repeat_unit": "ru-1", "extrapolated": {"x": 1}}
    )
    _seed_workflow(conn, "oligomer_sweep", "ru-2", status="error", summary={"repeat_unit": "ru-2"})
    _seed_workflow(conn, "oligomer_sweep", "ru-3")  # complete but NULL summary
    _seed_workflow(conn, "lambda_h", "ru-4", summary={"lambda_h_ev": 0.3})

    sweeps = db.list_complete_oligomer_sweeps(conn)
    assert len(sweeps) == 1
    assert sweeps[0]["workflow_id"] == ok_id
    assert sweeps[0]["smiles"] == "ru-1"
    assert sweeps[0]["extrapolated"] == {"x": 1}
    conn.close()


def test_candidate_summary_joins_three_sources(tmp_path):
    import json

    conn = db.open(str(tmp_path / "r.db"))
    smi = "c1ccccc1"
    jid = db.insert_job(conn, _job(signature="sig-p", smiles=smi, status="complete"))
    db.insert_result(conn, db.Result(job_id=jid, homo_ev=-5.2, lumo_ev=-2.9, gap_ev=2.3))

    tid = db.insert_job(
        conn, _job(signature="sig-t", smiles=smi, job_kind="tddft", status="complete")
    )
    states = [
        {"n": 1, "energy_ev": 2.5, "wavelength_nm": 495.9, "f": 0.10},
        {"n": 2, "energy_ev": 3.0, "wavelength_nm": 413.3, "f": 0.90},
        {"n": 3, "energy_ev": 3.5, "wavelength_nm": 354.2, "f": 0.30},
    ]
    db.insert_result(conn, db.Result(job_id=tid, spectra_json=json.dumps(states)))

    _seed_workflow(
        conn,
        "lambda_h",
        smi,
        summary={"lambda_h_ev": 0.31, "j_hole_ev": 0.05, "marcus_rate_s1": 1.2e12},
    )

    out = db.candidate_summary(conn, smi)
    assert out["homo_ev"] == -5.2 and out["lumo_ev"] == -2.9 and out["gap_ev"] == 2.3
    # Brightest state (middle one, f=0.90) wins, not the lowest.
    assert out["lambda_max_nm"] == 413.3 and out["f_osc"] == 0.90
    assert out["lambda_h_ev"] == 0.31
    assert out["j_hole_ev"] == 0.05
    assert out["marcus_rate_s1"] == 1.2e12
    conn.close()


def test_candidate_summary_missing_pieces_are_none(tmp_path):
    conn = db.open(str(tmp_path / "r.db"))
    smi = "c1ccccc1"
    jid = db.insert_job(conn, _job(signature="sig-p", smiles=smi, status="complete"))
    db.insert_result(conn, db.Result(job_id=jid, homo_ev=-5.2, lumo_ev=-2.9, gap_ev=2.3))

    out = db.candidate_summary(conn, smi)
    assert out["homo_ev"] == -5.2
    assert out["lambda_max_nm"] is None and out["f_osc"] is None
    assert out["lambda_h_ev"] is None and out["marcus_rate_s1"] is None
    conn.close()


def test_candidate_summary_ignores_cation_and_stale_rows(tmp_path):
    conn = db.open(str(tmp_path / "r.db"))
    smi = "c1ccccc1"
    old = db.insert_job(conn, _job(signature="sig-old", smiles=smi, status="complete"))
    db.insert_result(conn, db.Result(job_id=old, homo_ev=-9.9, lumo_ev=-9.9, gap_ev=0.0))
    new = db.insert_job(conn, _job(signature="sig-new", smiles=smi, status="complete"))
    db.insert_result(conn, db.Result(job_id=new, homo_ev=-5.2, lumo_ev=-2.9, gap_ev=2.3))
    # cation_opt: job_kind='properties' but charge=+1 — must not win even though newest.
    cat = db.insert_job(
        conn, _job(signature="sig-cat", smiles=smi, charge=1, mult=2, status="complete")
    )
    db.insert_result(conn, db.Result(job_id=cat, homo_ev=-11.0, lumo_ev=-8.0, gap_ev=3.0))

    out = db.candidate_summary(conn, smi)
    assert out["homo_ev"] == -5.2  # latest charge-0 row, not the cation, not the stale one
    conn.close()


def test_list_candidate_smiles_union(tmp_path):
    import json

    conn = db.open(str(tmp_path / "r.db"))
    db.insert_job(conn, _job(signature="s-a", smiles="A", status="complete"))
    tid = db.insert_job(conn, _job(signature="s-b", smiles="B", job_kind="tddft", status="complete"))
    db.insert_result(conn, db.Result(job_id=tid, spectra_json=json.dumps([])))
    _seed_workflow(conn, "lambda_h", "C", summary={"lambda_h_ev": 0.3})
    db.insert_job(conn, _job(signature="s-d", smiles="D", status="pending"))  # incomplete: excluded

    assert db.list_candidate_smiles(conn) == ["A", "B", "C"]
    conn.close()
