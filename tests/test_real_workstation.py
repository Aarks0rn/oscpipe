"""Real workstation integration tests — require user@203.0.113.10 with g16.

Run with:  pytest -m real
Skip by default in CI / offline runs.

The `real_settings` fixture is defined in conftest.py and skips automatically
if no SSH key is found.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import pytest

from oscpipe.cli.main import run_preflight, run_submit
from oscpipe.dft import gaussian
from oscpipe.store import db

_FIXTURES = Path(__file__).parent / "fixtures"


@pytest.mark.real
def test_real_preflight_passes(real_settings, capsys):
    """SSH connect, which g16, qstat, scratch writable — all green."""
    from oscpipe.dispatch.ssh import SshBackend

    conn = db.open(real_settings.db_path)
    backend = SshBackend(real_settings)
    try:
        rc = run_preflight(argparse.Namespace(), real_settings, backend, conn)
    finally:
        backend.close()

    assert rc == 0
    out = capsys.readouterr().out
    assert "g16: [ok]" in out
    assert "qstat: [ok]" in out
    assert "scratch: [ok]" in out


@pytest.mark.real
def test_real_submit_h2_full_cycle(real_settings, tmp_path):
    """Submit H2 (HF/STO-3G) via qsub; poll until complete; parse HOMO/LUMO."""
    import ase
    from oscpipe.dispatch.ssh import SshBackend

    backend = SshBackend(real_settings)
    atoms = ase.Atoms("H2", positions=[(0, 0, 0), (0, 0, 0.741)])
    com_text = gaussian.write_com_properties(
        atoms,
        "hf",
        "sto-3g",
        0,
        1,
        real_settings.gaussian_nproc,
        real_settings.gaussian_mem,
        "h2_pytest",
        "h2_pytest.chk",
    )

    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        com_path = Path(tmp) / "h2_pytest.com"
        com_path.write_text(com_text)
        jobid = backend.submit(str(com_path), "h2_pytest")

    deadline = time.time() + 300
    status = "pending"
    while time.time() < deadline:
        status = backend.poll(jobid)
        if status in ("complete", "error"):
            break
        time.sleep(real_settings.poll_interval_seconds)

    assert status == "complete", f"job timed out or errored (status={status})"

    log_path = backend.fetch_log(jobid, "h2_pytest", str(tmp_path / "logs"))
    backend.close()

    props = gaussian.parse_properties(log_path)
    # opt run → equilibrium geometry → HF/STO-3G HOMO ≈ -16.06 eV
    assert abs(props.homo_ev - (-16.060)) < 0.1
    assert props.gap_ev > 0
    assert abs(props.dipole_debye - 0.0) < 0.01
