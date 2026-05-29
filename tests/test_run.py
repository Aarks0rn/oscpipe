"""Tests for the JobRunner.resolve seam (cache → submit → persist-by-kind).

These exercise behaviour introduced when the single-shot CLI path and the λ_h
workflow were unified onto `JobRunner.resolve`:
  - a signature cache hit short-circuits without rebuilding the .com,
  - `need_log` recomputes when a cached job's local .log is gone,
  - persist-by-kind is total: sp_neutral/sp_cation get a results row, sp_dimer
    does not (its J_hole lives in the workflow summary),
  - re-running `oscpipe lambda` reuses completed sub-jobs from the cache.
"""

from __future__ import annotations

import argparse
import io
import json

from _helpers import StubBackend, synth_log

from oscpipe.cli.main import run_fetch, run_lambda
from oscpipe.run import JobRunner
from oscpipe.settings import Settings
from oscpipe.store import db
from oscpipe.store.cache import signature


def _settings(tmp_path):
    return Settings(
        backend="local",
        db_path=str(tmp_path / "results.db"),
        gaussian_nproc=2,
        gaussian_mem="2GB",
        poll_interval_seconds=0,
    )


def _job(sig, **overrides):
    base = dict(
        id=None,
        signature=sig,
        smiles="c1ccccc1",
        method="b3lyp",
        basis="6-31g**",
        charge=0,
        mult=1,
        job_kind="properties",
        status="pending",
        submitted_at="2026-05-29T00:00:00",
    )
    base.update(overrides)
    return db.Job(**base)


# λ_h per-suffix energies (match test_cli_lambda so λ_h ≠ 0, marcus stays finite).
_LAMBDA_ENERGIES_HA = {
    "neutral_opt": -154.0,
    "cation_opt": -153.5,
    "neutral_at_cation_geom": -153.95,
    "cation_at_neutral_geom": -153.45,
}


def _log_for_lambda(label: str) -> str:
    for suffix, energy in _LAMBDA_ENERGIES_HA.items():
        if label.endswith(suffix):
            return synth_log(energy_ha=energy)
    return synth_log()


# ── resolve: cache hit ───────────────────────────────────────────────────────


def test_resolve_cache_hit_returns_result_without_rebuilding(tmp_path):
    s = _settings(tmp_path)
    conn = db.open(s.db_path)
    sig = signature("c1ccccc1", "b3lyp", "6-31g**", 0, 1)
    jid = db.insert_job(conn, _job(sig, status="complete", log_path="/seeded.log"))
    db.insert_result(
        conn, db.Result(job_id=jid, homo_ev=-6.0, lumo_ev=-1.0, gap_ev=5.0, energy_ev=-4000.0)
    )

    runner = JobRunner(conn, s, StubBackend(), io.StringIO())
    built: list[int] = []

    def _build():
        built.append(1)
        return "com"

    r = runner.resolve(_job(sig), "lbl", _build, wait=False)

    assert r.cached is True
    assert r.result.energy_ev == -4000.0
    assert built == []  # .com was not rebuilt
    assert conn.execute("SELECT count(*) FROM jobs").fetchone()[0] == 1  # no new submit


def test_resolve_need_log_recomputes_when_cached_log_missing(tmp_path):
    s = _settings(tmp_path)
    conn = db.open(s.db_path)
    sig = signature("c1ccccc1", "b3lyp", "6-31g**", 0, 1)
    jid = db.insert_job(conn, _job(sig, status="complete", log_path=str(tmp_path / "gone.log")))
    db.insert_result(
        conn, db.Result(job_id=jid, homo_ev=-6.0, lumo_ev=-1.0, gap_ev=5.0, energy_ev=-4000.0)
    )

    backend = StubBackend(log_text=synth_log(energy_ha=-100.0))
    backend.poll_return = "complete"
    runner = JobRunner(conn, s, backend, io.StringIO())
    built: list[int] = []

    def _build():
        built.append(1)
        return "com"

    r = runner.resolve(_job(sig), "lbl", _build, wait=True, need_log=True)

    assert built == [1]  # cached .log gone → recomputed
    assert r.cached is False
    assert r.status == "complete"
    assert conn.execute("SELECT count(*) FROM jobs").fetchone()[0] == 2


# ── persist-by-kind totality (through λ_h) ───────────────────────────────────


def test_lambda_persists_sp_results_but_not_dimer(tmp_path):
    s = _settings(tmp_path)
    backend = StubBackend(log_provider=_log_for_lambda)
    backend.poll_return = "complete"
    conn = db.open(s.db_path)
    from oscpipe.chem.smiles import embed_3d

    initial = embed_3d("c1ccccc1")
    rc = run_lambda(
        argparse.Namespace(smiles="c1ccccc1"), s, backend, conn, geometry_loader=lambda _: initial
    )
    assert rc == 0

    for kind in ("sp_neutral", "sp_cation"):
        jrow = conn.execute("SELECT id FROM jobs WHERE job_kind = ?", (kind,)).fetchone()
        res = conn.execute(
            "SELECT energy_ev FROM results WHERE job_id = ?", (jrow["id"],)
        ).fetchone()
        assert res is not None, f"{kind} should have a results row"
        assert res["energy_ev"] is not None

    drow = conn.execute("SELECT id FROM jobs WHERE job_kind = 'sp_dimer'").fetchone()
    assert conn.execute("SELECT 1 FROM results WHERE job_id = ?", (drow["id"],)).fetchone() is None


def test_lambda_rerun_reuses_cache_no_new_jobs(tmp_path):
    s = _settings(tmp_path)
    backend = StubBackend(log_provider=_log_for_lambda)
    backend.poll_return = "complete"
    conn = db.open(s.db_path)
    from oscpipe.chem.smiles import embed_3d

    initial = embed_3d("c1ccccc1")
    args = argparse.Namespace(smiles="c1ccccc1")

    assert run_lambda(args, s, backend, conn, geometry_loader=lambda _: initial) == 0
    n_after_first = conn.execute("SELECT count(*) FROM jobs").fetchone()[0]
    assert n_after_first == 5

    # Second run: every sub-job signature hits the cache → no new jobs submitted.
    assert run_lambda(args, s, backend, conn, geometry_loader=lambda _: initial) == 0
    assert conn.execute("SELECT count(*) FROM jobs").fetchone()[0] == 5
    # Two workflows recorded, both complete.
    wfs = list(
        conn.execute(
            "SELECT status, summary_json FROM workflows WHERE kind = 'lambda_h' ORDER BY id"
        )
    )
    assert len(wfs) == 2
    assert all(w["status"] == "complete" for w in wfs)
    # The cached re-run reproduces the same λ_h, not just the same job count.
    second = json.loads(wfs[1]["summary_json"])
    assert abs(second["lambda_h_ev"] - 2.72114) < 1e-2
    assert second["j_hole_ev"] is not None


# ── persist-by-kind on the resume (fetch) path ───────────────────────────────


def test_fetch_resumes_dimer_without_bogus_results_row(tmp_path):
    s = _settings(tmp_path)
    conn = db.open(s.db_path)
    sig = signature("c1ccccc1", "wb97xd", "6-31g**", 0, 1, extras="dimer_sp")
    jid = db.insert_job(
        conn, _job(sig, method="wb97xd", job_kind="sp_dimer", status="running", ssh_jobid="d-1")
    )

    backend = StubBackend()  # default SYNTH_LOG has orbitals — would mis-persist as properties
    backend.poll_return = "complete"
    run_fetch(argparse.Namespace(job_id=None), s, backend, conn)

    row = conn.execute("SELECT status FROM jobs WHERE id = ?", (jid,)).fetchone()
    assert row["status"] == "complete"
    assert conn.execute("SELECT 1 FROM results WHERE job_id = ?", (jid,)).fetchone() is None


# ── Properties browser filter (mirrors app/pages/1_Properties.py) ─────────────


def test_properties_browser_hides_lambda_subjobs(tmp_path):
    """Total-persist gives λ_h opt jobs (job_kind='properties') a results row; the
    Properties browser must hide those workflow internals while still showing
    standalone and batch property jobs."""
    s = _settings(tmp_path)
    conn = db.open(s.db_path)
    lam_wf = db.insert_workflow(conn, "lambda_h", "c1ccccc1", "2026-05-29T00:00:00")
    batch_wf = db.insert_workflow(conn, "batch", "batch:1", "2026-05-29T00:00:00")

    def _complete_prop(smiles, sig, workflow_id):
        jid = db.insert_job(
            conn, _job(sig, smiles=smiles, status="complete", workflow_id=workflow_id)
        )
        db.insert_result(
            conn,
            db.Result(job_id=jid, homo_ev=-6.0, lumo_ev=-1.0, gap_ev=5.0, energy_ev=-4000.0),
        )
        return jid

    lam_job = _complete_prop("c1ccccc1", "sig-lam", lam_wf)
    standalone = _complete_prop("c1ccsc1", "sig-std", None)
    batch_job = _complete_prop("c1ccoc1", "sig-bat", batch_wf)

    rows = conn.execute(
        "SELECT id FROM v_jobs_with_results WHERE job_kind = 'properties' "
        "AND (workflow_id IS NULL OR workflow_id NOT IN "
        "(SELECT id FROM workflows WHERE kind = 'lambda_h')) "
        "ORDER BY id"
    ).fetchall()
    ids = {r["id"] for r in rows}
    assert standalone in ids  # standalone submit shown
    assert batch_job in ids  # batch job (has workflow_id, but kind=batch) shown
    assert lam_job not in ids  # λ_h sub-job hidden from the properties browser
