"""Tests for `oscpipe oligomer` — oligomer-length sweep + 1/n extrapolation.

Per length n: a B3LYP `properties` (opt) job for HOMO/gap and an ωB97X-D `tddft`
job for the optical gap, then each property is extrapolated vs 1/n to the polymer
limit. The StubBackend returns per-n synthetic logs so the extrapolated limits
are predictable: HOMO_ev = -4.0 - 1.0*(1/n), optical_gap_ev = 1.8 + 1.2*(1/n),
so the limits at 1/n -> 0 are -4.0 and 1.8 eV.
"""

from __future__ import annotations

import argparse
import json
import re

import ase
import pytest
from _helpers import StubBackend

from oscpipe.cli.main import run_oligomer
from oscpipe.settings import Settings
from oscpipe.store import db

HARTREE_TO_EV = 27.2114


def _props_log(homo_ev: float, lumo_ev: float = -2.0) -> str:
    occ = [(homo_ev - 0.3) / HARTREE_TO_EV, homo_ev / HARTREE_TO_EV]
    virt = [lumo_ev / HARTREE_TO_EV, (lumo_ev + 0.5) / HARTREE_TO_EV]
    return (
        " SCF Done:  E(RB3LYP) =       -154.000000 A.U. after   8 cycles\n"
        f" Alpha  occ. eigenvalues --   {occ[0]:.5f}  {occ[1]:.5f}\n"
        f" Alpha virt. eigenvalues --    {virt[0]:.5f}   {virt[1]:.5f}\n"
        " Dipole moment (field-independent basis, Debye):\n"
        "    X=    0.0000  Y=    0.0000  Z=    0.0000  Tot=    0.0000\n"
        " Normal termination of Gaussian 16 at ...\n"
    )


def _tddft_log(s1_ev: float) -> str:
    return (
        " ## stub TDDFT\n"
        f" Excited State   1:      Singlet-A   {s1_ev:.4f} eV  {1239.8 / s1_ev:.2f} nm  f=1.5000\n"
        f" Excited State   2:      Singlet-A   {s1_ev + 0.6:.4f} eV  "
        f"{1239.8 / (s1_ev + 0.6):.2f} nm  f=0.3000\n"
        " Normal termination of Gaussian 16 at ...\n"
    )


def _sweep_log_provider(label: str) -> str:
    m = re.search(r"n(\d+)_(opt|tddft)$", label)
    assert m, f"unexpected label: {label}"
    n, kind = int(m.group(1)), m.group(2)
    if kind == "opt":
        return _props_log(homo_ev=-4.0 - 1.0 / n)
    return _tddft_log(s1_ev=1.8 + 1.2 / n)


def _args(repeat_unit="[*]c1ccc([*])s1", max_n=3, nstates=5):
    return argparse.Namespace(repeat_unit=repeat_unit, max_n=max_n, nstates=nstates)


def _setup(tmp_path):
    s = Settings(
        backend="local",
        db_path=str(tmp_path / "results.db"),
        gaussian_nproc=2,
        gaussian_mem="2GB",
        poll_interval_seconds=0,
    )
    backend = StubBackend(log_provider=_sweep_log_provider)
    backend.poll_return = "complete"
    conn = db.open(s.db_path)
    return s, backend, conn


def test_sweep_runs_properties_and_tddft_for_each_length(tmp_path):
    s, backend, conn = _setup(tmp_path)
    stub_geom = ase.Atoms("C", positions=[[0.0, 0.0, 0.0]])

    rc = run_oligomer(_args(max_n=3), s, backend, conn, geometry_loader=lambda _: stub_geom)
    assert rc == 0

    jobs = list(conn.execute("SELECT * FROM jobs ORDER BY id"))
    assert len(jobs) == 6  # 3 lengths x (properties + tddft)
    assert all(j["status"] == "complete" for j in jobs)
    assert len({j["workflow_id"] for j in jobs}) == 1
    kinds = sorted(j["job_kind"] for j in jobs)
    assert kinds == ["properties", "properties", "properties", "tddft", "tddft", "tddft"]


def test_sweep_extrapolates_to_polymer_limit(tmp_path):
    s, backend, conn = _setup(tmp_path)
    stub_geom = ase.Atoms("C", positions=[[0.0, 0.0, 0.0]])

    rc = run_oligomer(_args(max_n=3), s, backend, conn, geometry_loader=lambda _: stub_geom)
    assert rc == 0

    wf = conn.execute("SELECT * FROM workflows").fetchone()
    assert wf["status"] == "complete"
    summary = json.loads(wf["summary_json"])
    assert len(summary["per_n"]) == 3
    # Limits recovered to within the 5-decimal round-trip noise of the synthetic logs.
    assert summary["extrapolated"]["homo_ev"]["limit"] == pytest.approx(-4.0, abs=1e-2)
    assert summary["extrapolated"]["optical_gap_ev"]["limit"] == pytest.approx(1.8, abs=1e-2)
    assert summary["extrapolated"]["optical_gap_ev"]["r_squared"] == pytest.approx(1.0, abs=1e-4)


def test_sweep_opt_failure_marks_workflow_error(tmp_path):
    s, backend, conn = _setup(tmp_path)
    backend.poll_return = "error"

    rc = run_oligomer(_args(max_n=2), s, backend, conn, geometry_loader=lambda _: None)
    assert rc == 1
    wf = conn.execute("SELECT * FROM workflows").fetchone()
    assert wf["status"] == "error"
