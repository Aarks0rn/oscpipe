"""Tests for `oscpipe lambda` (4-point λ_h workflow)."""

from __future__ import annotations

import argparse
import json

from _helpers import StubBackend, synth_log

from oscpipe.cli.main import run_lambda
from oscpipe.settings import Settings
from oscpipe.store import db

# Energies in Hartree assigned per label suffix. Convergent stabilisation:
#   neutral_opt at neutral geom         = -154.0
#   cation_opt  at cation geom          = -153.5
#   neutral SP  at cation geom (worse)  = -153.95
#   cation  SP  at neutral geom (worse) = -153.45
#
# λ_h = (E_cat(N) − E_cat(C)) + (E_neut(C) − E_neut(N))
#     = (−153.45 − (−153.5)) + (−153.95 − (−154.0)) = 0.10 Ha
#     ≈ 2.72114 eV
ENERGIES_HA = {
    "neutral_opt": -154.0,
    "cation_opt": -153.5,
    "neutral_at_cation_geom": -153.95,
    "cation_at_neutral_geom": -153.45,
}


def _log_for(label: str) -> str:
    for suffix, energy in ENERGIES_HA.items():
        if label.endswith(suffix):
            return synth_log(energy_ha=energy)
    return synth_log()  # default


def _args(smiles: str = "c1ccccc1"):
    return argparse.Namespace(smiles=smiles)


def _setup(tmp_path):
    s = Settings(
        backend="local",
        db_path=str(tmp_path / "results.db"),
        gaussian_nproc=2,
        gaussian_mem="2GB",
        poll_interval_seconds=0,
    )
    backend = StubBackend(log_provider=_log_for)
    backend.poll_return = "complete"
    conn = db.open(s.db_path)
    return s, backend, conn


def test_lambda_inserts_five_jobs_with_same_workflow_id(tmp_path):
    s, backend, conn = _setup(tmp_path)
    initial = None

    def loader(_log_path):
        # geometry loader called after each opt; reuse the initial embed.
        nonlocal initial
        if initial is None:
            from oscpipe.chem.smiles import embed_3d

            initial = embed_3d("c1ccccc1")
        return initial

    rc = run_lambda(_args(), s, backend, conn, geometry_loader=loader)
    assert rc == 0

    jobs = list(conn.execute("SELECT * FROM jobs ORDER BY id"))
    assert len(jobs) == 5
    wf_ids = {j["workflow_id"] for j in jobs}
    assert len(wf_ids) == 1
    assert all(j["status"] == "complete" for j in jobs)
    kinds = {j["job_kind"] for j in jobs}
    assert kinds == {"properties", "sp_neutral", "sp_cation", "sp_dimer"}


def test_lambda_computes_lambda_h(tmp_path, capsys):
    s, backend, conn = _setup(tmp_path)
    from oscpipe.chem.smiles import embed_3d

    initial = embed_3d("c1ccccc1")
    rc = run_lambda(_args(), s, backend, conn, geometry_loader=lambda _: initial)
    assert rc == 0
    wf = conn.execute("SELECT * FROM workflows").fetchone()
    assert wf["status"] == "complete"
    summary = json.loads(wf["summary_json"])
    # λ_h = 0.10 Ha × 27.2114 = 2.72114 eV
    assert abs(summary["lambda_h_ev"] - 2.72114) < 1e-2
    # Sanity-check stored energies.
    e = summary["energies_ev"]
    HARTREE_TO_EV = 27.2114
    assert abs(e["neutral_opt"] - (-154.0 * HARTREE_TO_EV)) < 1e-2
    assert abs(e["cation_opt"] - (-153.5 * HARTREE_TO_EV)) < 1e-2
    # Dimer SP succeeded → J_hole and marcus_rate_s1 present.
    assert summary["j_hole_ev"] is not None
    assert summary["j_hole_ev"] >= 0.0
    assert summary["marcus_rate_s1"] is not None
    assert summary["marcus_rate_s1"] > 0.0
    out = capsys.readouterr().out
    assert "lambda_h = 2.7211 eV" in out
    assert "J_hole" in out
    assert "marcus_rate" in out


def test_lambda_opt_stage_failure_marks_workflow_error(tmp_path, capsys):
    s, backend, conn = _setup(tmp_path)
    backend.poll_return = "error"

    rc = run_lambda(_args(), s, backend, conn, geometry_loader=lambda _: None)
    assert rc == 1
    wf = conn.execute("SELECT * FROM workflows").fetchone()
    assert wf["status"] == "error"
    summary = json.loads(wf["summary_json"])
    assert summary["stage"] == "opt"
    assert "lambda: opt stage failed" in capsys.readouterr().out
