"""Tests for `oscpipe uvvis` (TDDFT submit + excited-state parse)."""

from __future__ import annotations

import argparse
import json

from _helpers import StubBackend, synth_tddft_log

from oscpipe.cli.main import run_uvvis
from oscpipe.settings import Settings
from oscpipe.store import db


def _args(smiles: str, nstates: int = 3, method: str = "b3lyp", from_log=None):
    return argparse.Namespace(
        smiles=smiles, nstates=nstates, method=method, from_log=from_log
    )


def _setup(tmp_path):
    s = Settings(
        backend="local",
        db_path=str(tmp_path / "results.db"),
        gaussian_nproc=2,
        gaussian_mem="2GB",
    )
    backend = StubBackend(log_text=synth_tddft_log(n=3))
    backend.poll_return = "complete"
    conn = db.open(s.db_path)
    return s, backend, conn


def test_uvvis_inserts_tddft_job_and_parses_states(tmp_path, capsys):
    s, backend, conn = _setup(tmp_path)
    rc = run_uvvis(_args("c1ccccc1", nstates=3), s, backend, conn)
    assert rc == 0

    row = conn.execute("SELECT * FROM jobs ORDER BY id DESC LIMIT 1").fetchone()
    assert row["job_kind"] == "tddft"
    assert row["status"] == "complete"
    assert row["notes"] == "nstates=3"

    result = conn.execute(
        "SELECT spectra_json FROM results WHERE job_id = ?", (row["id"],)
    ).fetchone()
    states = json.loads(result["spectra_json"])
    assert len(states) == 3
    assert {"n", "energy_ev", "wavelength_nm", "f"} <= states[0].keys()
    assert abs(states[0]["energy_ev"] - 2.5) < 1e-3
    assert abs(states[2]["energy_ev"] - 3.5) < 1e-3

    out = capsys.readouterr().out
    assert "complete: job=" in out
    assert "tddft n=3" in out


def test_uvvis_cache_hit_short_circuits(tmp_path, capsys):
    s, backend, conn = _setup(tmp_path)
    # First call populates cache.
    run_uvvis(_args("c1ccccc1", nstates=3), s, backend, conn)
    n_jobs_before = conn.execute("SELECT count(*) FROM jobs").fetchone()[0]

    # Second call should hit cache.
    capsys.readouterr()
    rc = run_uvvis(_args("c1ccccc1", nstates=3), s, backend, conn)
    assert rc == 0
    n_jobs_after = conn.execute("SELECT count(*) FROM jobs").fetchone()[0]
    assert n_jobs_after == n_jobs_before
    assert "cache hit" in capsys.readouterr().out


def test_uvvis_different_nstates_no_cache_collide(tmp_path):
    s, backend, conn = _setup(tmp_path)
    run_uvvis(_args("c1ccccc1", nstates=3), s, backend, conn)
    # Same SMILES, different nstates → different signature → new submit.
    run_uvvis(_args("c1ccccc1", nstates=5), s, backend, conn)
    n = conn.execute("SELECT count(*) FROM jobs").fetchone()[0]
    assert n == 2


def test_uvvis_from_log_reuses_geometry_skipping_embed(tmp_path, capsys):
    """--from-log path: geometry_loader fires, RDKit embed is bypassed."""
    s, backend, conn = _setup(tmp_path)
    fake_log = tmp_path / "neutral_opt.log"
    fake_log.write_text("stub")

    calls: list[tuple[str, str]] = []

    def loader(log_path: str, canonical: str):
        calls.append((log_path, canonical))
        from oscpipe.chem.smiles import embed_3d

        return embed_3d(canonical)

    rc = run_uvvis(
        _args("c1ccccc1", nstates=3, from_log=str(fake_log)),
        s,
        backend,
        conn,
        geometry_loader=loader,
    )
    assert rc == 0
    assert len(calls) == 1 and calls[0][0] == str(fake_log)
    out = capsys.readouterr().out
    assert "reusing geometry from" in out


def test_uvvis_from_log_validator_rejects_atom_mismatch(tmp_path):
    """_atoms_from_log_validated raises on element-multiset mismatch."""
    from oscpipe.cli.main import _atoms_from_log_validated

    # Build a "log" whose ASE-read result is methane (CH4); SMILES says benzene.
    fake_log = tmp_path / "wrong_molecule.log"
    fake_log.write_text("stub")

    import pytest
    from oscpipe.chem.smiles import embed_3d

    def fake_reader(_path):
        return embed_3d("C")  # methane: 1 C, 4 H

    import oscpipe.cli.main as cli_main

    real_reader = cli_main.read_gaussian_log
    cli_main.read_gaussian_log = fake_reader
    try:
        with pytest.raises(ValueError, match="atom multiset mismatch"):
            _atoms_from_log_validated(str(fake_log), "c1ccccc1")  # benzene: 6 C, 6 H
    finally:
        cli_main.read_gaussian_log = real_reader


def test_uvvis_from_log_validator_accepts_matching_atoms(tmp_path):
    """_atoms_from_log_validated passes when atom multiset matches the SMILES."""
    from oscpipe.cli.main import _atoms_from_log_validated

    fake_log = tmp_path / "benzene_opt.log"
    fake_log.write_text("stub")

    from oscpipe.chem.smiles import embed_3d

    def fake_reader(_path):
        return embed_3d("c1ccccc1")

    import oscpipe.cli.main as cli_main

    real_reader = cli_main.read_gaussian_log
    cli_main.read_gaussian_log = fake_reader
    try:
        atoms = _atoms_from_log_validated(str(fake_log), "c1ccccc1")
        assert atoms is not None
        assert len(atoms.get_chemical_symbols()) > 0
    finally:
        cli_main.read_gaussian_log = real_reader


def test_uvvis_from_log_missing_file_raises(tmp_path):
    """_atoms_from_log_validated raises FileNotFoundError for a missing log."""
    from oscpipe.cli.main import _atoms_from_log_validated

    import pytest

    with pytest.raises(FileNotFoundError, match="does not exist"):
        _atoms_from_log_validated(str(tmp_path / "nope.log"), "c1ccccc1")


def test_uvvis_does_not_collide_with_properties_cache(tmp_path):
    """A pre-existing properties job must not satisfy a uvvis cache lookup."""
    s, backend, conn = _setup(tmp_path)
    # Seed a complete properties job with the same SMILES / method / basis.
    from oscpipe.chem.smiles import canonicalise
    from oscpipe.store.cache import signature

    canon, _ = canonicalise("c1ccccc1")
    prop_sig = signature(canon, "b3lyp", "6-31g*", 0, 1)  # job_kind=properties
    jid = db.insert_job(
        conn,
        db.Job(
            id=None,
            signature=prop_sig,
            smiles=canon,
            method="b3lyp",
            basis="6-31g*",
            charge=0,
            mult=1,
            job_kind="properties",
            status="complete",
            submitted_at="2026-05-21T00:00:00",
        ),
    )
    db.insert_result(conn, db.Result(job_id=jid, homo_ev=-6.7, lumo_ev=-1.0, gap_ev=5.7))

    run_uvvis(_args("c1ccccc1", nstates=3), s, backend, conn)
    # uvvis should have submitted a separate tddft row.
    n_tddft = conn.execute("SELECT count(*) FROM jobs WHERE job_kind='tddft'").fetchone()[0]
    assert n_tddft == 1
