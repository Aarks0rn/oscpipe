"""Tests for `oscpipe lambda` (4-point λ_h workflow)."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import ase
import numpy as np
from _helpers import StubBackend, synth_log

from oscpipe.chem.geometry import build_pi_stack_dimer
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


class _ConcBackend:
    """Tracks peak concurrent in-flight jobs. Each job polls 'running' once, then
    'complete', so jobs in the same layer genuinely overlap."""

    def __init__(self, log_provider):
        self.active = 0
        self.max_active = 0
        self._polls: dict[str, int] = {}
        self._label: dict[str, str] = {}
        self._jobs: dict[str, str] = {}
        self.log_provider = log_provider

    def submit(self, com_path, label):
        sid = f"stub-{label}"
        self.active += 1
        self.max_active = max(self.max_active, self.active)
        self._polls[sid] = 0
        self._label[sid] = label
        return sid

    def poll(self, sid):
        self._polls[sid] += 1
        return "running" if self._polls[sid] <= 1 else "complete"

    def fetch_log(self, sid, label, local_dir):
        self.active -= 1
        Path(local_dir).mkdir(parents=True, exist_ok=True)
        out = Path(local_dir) / f"{label}.log"
        out.write_text(self.log_provider(label))
        return str(out)

    def cancel(self, sid):
        pass


def test_lambda_runs_opts_in_parallel_and_matches_serial(tmp_path):
    s = Settings(
        backend="local",
        db_path=str(tmp_path / "results.db"),
        gaussian_nproc=10,
        gaussian_mem="24GB",
        poll_interval_seconds=0,
        max_lanes=2,
    )
    backend = _ConcBackend(_log_for)
    conn = db.open(s.db_path)
    from oscpipe.chem.smiles import embed_3d

    initial = embed_3d("c1ccccc1")
    rc = run_lambda(_args(), s, backend, conn, geometry_loader=lambda _: initial)
    assert rc == 0
    # Two lanes were actually used — the neutral/cation opts overlapped.
    assert backend.max_active == 2
    summary = json.loads(
        conn.execute("SELECT summary_json FROM workflows").fetchone()["summary_json"]
    )
    # Same physics as the serial path (test_lambda_computes_lambda_h).
    assert abs(summary["lambda_h_ev"] - 2.72114) < 1e-2
    assert summary["j_hole_ev"] is not None


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


def test_build_pi_stack_dimer_no_clash_when_misaligned():
    """Stacked dimer must have no clashing atoms even if the monomer is not XY-aligned."""
    # Four coplanar C atoms in the XY plane, then rotated 60° around X.
    # Old code: second copy at Z_original + 3.5 → overlap. New code: aligns first.
    atoms = ase.Atoms(
        "CCCC",
        positions=[[0.0, 0.0, 0.0], [1.4, 0.0, 0.0], [2.8, 0.0, 0.0], [4.2, 0.0, 0.0]],
    )
    theta = np.radians(60)
    Rx = np.array(
        [[1, 0, 0], [0, np.cos(theta), -np.sin(theta)], [0, np.sin(theta), np.cos(theta)]]
    )
    atoms.positions = (Rx @ atoms.positions.T).T

    stack_distance = 3.5
    dimer = build_pi_stack_dimer(atoms, stack_distance=stack_distance)

    n = len(atoms)
    pos1 = dimer.positions[:n]
    pos2 = dimer.positions[n:]

    dists = np.linalg.norm(pos1[:, np.newaxis] - pos2[np.newaxis], axis=2)
    assert dists.min() > 2.5, f"atoms too close: {dists.min():.2f} Å"
    assert abs(pos2[:, 2].mean() - pos1[:, 2].mean() - stack_distance) < 0.05


def test_build_pi_stack_dimer_no_clash_for_nonplanar_molecule():
    """Dimer must not clash when the monomer has genuine out-of-plane atoms.

    Reproduces the real failure mode for C40H18F4N4S8: ring-twisted molecules
    have H atoms ±2.2 Å from the mean plane.  Old code (simple Z shift) puts
    the upper H of monomer1 at +2.2 Å and the lower H of monomer2 at 3.5-2.2=1.3 Å
    — only 0.9 Å apart.  New code must slip-stack until gap >= 2.5 Å.
    """
    # Four aromatic-core atoms near-planar, plus one atom far above and one far below.
    # This is the smallest geometry that reproduces the ±2.2 Å Z-extent of the real molecule.
    positions = [
        [0.0, 0.0, 0.0],
        [1.4, 0.0, 0.0],
        [0.0, 1.4, 0.0],
        [1.4, 1.4, 0.0],
        [0.7, 0.7, 2.2],  # far-above atom (H-like)
        [0.7, 0.7, -2.2],  # far-below atom (H-like)
    ]
    atoms = ase.Atoms("CCCCCC", positions=positions)

    stack_distance = 3.5
    dimer = build_pi_stack_dimer(atoms, stack_distance=stack_distance)

    n = len(atoms)
    pos1 = dimer.positions[:n]
    pos2 = dimer.positions[n:]

    dists = np.linalg.norm(pos1[:, np.newaxis] - pos2[np.newaxis], axis=2)
    assert dists.min() > 2.5, f"atoms too close: {dists.min():.2f} Å"
