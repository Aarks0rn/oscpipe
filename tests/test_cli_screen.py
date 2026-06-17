"""Tests for `oscpipe screen` (campaign: properties + tddft + λ_h per candidate)."""

from __future__ import annotations

import argparse
import io
import json

from _helpers import StubBackend, synth_log, synth_tddft_log

from oscpipe.cli.main import run_screen
from oscpipe.settings import Settings
from oscpipe.store import db

# λ_h sub-job energies per label suffix (same scheme as test_cli_lambda):
# λ_h = 0.10 Ha ≈ 2.72114 eV.  neutral_opt is NOT in this map — screen's own
# properties job doubles as λ_h's neutral_opt via the signature cache, so its
# energy comes from the combined log below.
ENERGIES_HA = {
    "cation_opt": -153.5,
    "neutral_at_cation_geom": -153.95,
    "cation_at_neutral_geom": -153.45,
}


def _log_for(label: str) -> str:
    for suffix, energy in ENERGIES_HA.items():
        if label.endswith(suffix):
            return synth_log(energy_ha=energy)
    # Screen's properties and tddft labels carry no suffix, so one combined log
    # must satisfy both parsers (and the dimer SP's eigenvalue read).  Its SCF
    # energy is neutral_opt's −154.0 Ha because λ_h reads e_neut_opt from this
    # job through the cache.
    return synth_log(energy_ha=-154.0) + synth_tddft_log(n=3)


def _args(file: str) -> argparse.Namespace:
    return argparse.Namespace(file=file, method="b3lyp", basis="6-31g**", nstates=3)


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


def _loader():
    """Deterministic geometry loader for the λ_h leg."""
    from oscpipe.chem.smiles import embed_3d

    initial = embed_3d("c1ccccc1")
    return lambda _log_path: initial


def _csv(tmp_path, smiles_rows):
    p = tmp_path / "candidates.csv"
    p.write_text("smiles\n" + "\n".join(smiles_rows) + "\n")
    return str(p)


def _run(args, s, backend, conn, stdin=None):
    # screen's workflow needs the injected loader; go through the workflow
    # directly for that, via the cli adapter signature.
    from oscpipe.cli.main import _read_smiles
    from oscpipe.workflows.screen import run_screen as run_screen_workflow

    smiles_list = _read_smiles(args.file, stdin or io.StringIO(""))
    return run_screen_workflow(
        smiles_list, args, s, backend, conn, stdout=io.StringIO(), geometry_loader=_loader()
    )


def test_screen_runs_three_stages_per_candidate(tmp_path):
    s, backend, conn = _setup(tmp_path)
    rc = _run(_args(_csv(tmp_path, ["c1ccccc1"])), s, backend, conn)
    assert rc == 0

    wf = conn.execute("SELECT * FROM workflows WHERE kind='screen'").fetchone()
    assert wf["status"] == "complete"
    summary = json.loads(wf["summary_json"])
    assert summary["total"] == 1 and summary["ok"] == 1 and summary["failed"] == 0
    cand = summary["candidates"][0]
    assert cand["properties"] == "ok"
    assert cand["tddft"] == "ok"
    assert cand["lambda_h"]["status"] == "complete"
    assert isinstance(cand["lambda_h"]["workflow_id"], int)

    # 6 jobs total: properties + tddft from screen, then cation_opt + 2 SPs +
    # dimer from λ_h — its neutral_opt is a cache hit of screen's properties job.
    n_jobs = conn.execute("SELECT count(*) FROM jobs").fetchone()[0]
    assert n_jobs == 6

    # λ_h produced the known value from the shared energies.
    lam = conn.execute("SELECT * FROM workflows WHERE kind='lambda_h'").fetchone()
    lam_summary = json.loads(lam["summary_json"])
    assert abs(lam_summary["lambda_h_ev"] - 2.72114) < 1e-2


def test_screen_continues_after_candidate_failure(tmp_path):
    s, backend, conn = _setup(tmp_path)
    # "!!" survives canonicalise (warning only) and then fails in embed_3d.
    rc = _run(_args(_csv(tmp_path, ["!!", "c1ccccc1"])), s, backend, conn)
    assert rc == 1

    wf = conn.execute("SELECT * FROM workflows WHERE kind='screen'").fetchone()
    assert wf["status"] == "error"
    summary = json.loads(wf["summary_json"])
    assert summary["failed"] == 1 and summary["ok"] == 1
    assert "error" in summary["candidates"][0]
    assert summary["candidates"][1]["lambda_h"]["status"] == "complete"


def test_screen_rerun_reuses_cache_no_new_jobs(tmp_path):
    s, backend, conn = _setup(tmp_path)
    args = _args(_csv(tmp_path, ["c1ccccc1"]))
    assert _run(args, s, backend, conn) == 0
    n_before = conn.execute("SELECT count(*) FROM jobs").fetchone()[0]

    assert _run(args, s, backend, conn) == 0
    n_after = conn.execute("SELECT count(*) FROM jobs").fetchone()[0]
    assert n_after == n_before

    second = conn.execute(
        "SELECT * FROM workflows WHERE kind='screen' ORDER BY id DESC LIMIT 1"
    ).fetchone()
    summary = json.loads(second["summary_json"])
    cand = summary["candidates"][0]
    assert cand["properties"] == "cached"
    assert cand["tddft"] == "cached"
    assert cand["lambda_h"]["status"] == "complete"


def test_screen_adapter_reads_csv_and_stdin(tmp_path):
    """The cli adapter feeds both file and stdin input into the workflow."""
    s, backend, conn = _setup(tmp_path)
    out = io.StringIO()
    rc = run_screen(
        _args("-"),
        s,
        backend,
        conn,
        stdin=io.StringIO("c1ccccc1\n"),
        stdout=out,
    )
    # No injected geometry_loader here: λ_h reads the real (combined) logs, so
    # only the adapter plumbing is being asserted, not the physics.
    assert "screen: workflow_id=" in out.getvalue()
    assert rc in (0, 1)

    out = io.StringIO()
    rc = run_screen(_args(_csv(tmp_path, [])), s, backend, conn, stdout=out)
    assert rc == 1
    assert "no SMILES read" in out.getvalue()
