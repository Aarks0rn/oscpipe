"""Shared test helpers — StubBackend, synthetic Gaussian log, DB seeders."""

from __future__ import annotations

from pathlib import Path

from oscpipe.store import db

SYNTH_LOG = """\
 SCF Done:  E(RB3LYP) =       -1.17000000000 A.U. after   8 cycles
 Alpha  occ. eigenvalues --   -0.55000  -0.50000
 Alpha virt. eigenvalues --    0.10000   0.20000
 Dipole moment (field-independent basis, Debye):
    X=    0.0000  Y=    0.0000  Z=    0.0000  Tot=    0.0000
 Normal termination of Gaussian 16 at Mon Jan  1 00:00:00 2026.
"""

EMPTY_LOG = " Normal termination of Gaussian 16 at ...\n"


def synth_log(energy_ha: float = -1.17) -> str:
    """Minimal Gaussian-format log with a parseable SCF energy."""
    return (
        f" SCF Done:  E(RB3LYP) =       {energy_ha:.6f} A.U. after   8 cycles\n"
        " Alpha  occ. eigenvalues --   -0.55000  -0.50000\n"
        " Alpha virt. eigenvalues --    0.10000   0.20000\n"
        " Dipole moment (field-independent basis, Debye):\n"
        "    X=    0.0000  Y=    0.0000  Z=    0.0000  Tot=    0.0000\n"
        " Normal termination of Gaussian 16 at ...\n"
    )


def synth_tddft_log(n: int = 3) -> str:
    """Minimal TDDFT log with `n` Excited State lines."""
    rows = []
    for i in range(1, n + 1):
        energy = 2.0 + 0.5 * i
        wavelength = 1239.8 / energy
        f = 0.1 * i
        rows.append(
            f" Excited State{i:>4}:      Singlet-A       {energy:.4f} eV "
            f"{wavelength:.2f} nm  f={f:.4f}  <S**2>=0.000"
        )
    return (
        " ## stub TDDFT log\n" + "\n".join(rows) + "\n Normal termination of Gaussian 16 at ...\n"
    )


class StubBackend:
    """Configurable Backend stand-in. Tests set the return values they need.

    Set `log_provider` to a callable(label) -> str for per-label logs (used
    by the λ_h workflow test, which needs distinct energies per sub-job).
    """

    def __init__(self, log_text: str = SYNTH_LOG, log_provider=None):
        self.poll_return = "running"
        self.log_text = log_text
        self.log_provider = log_provider
        self.calls: list[tuple[str, str]] = []
        self._jobs: dict[str, str] = {}

    def submit(self, com_path, label):
        self.calls.append(("submit", label))
        return f"stub-{label[:20]}"

    def poll(self, remote_job_id):
        self.calls.append(("poll", remote_job_id))
        return self.poll_return

    def fetch_log(self, remote_job_id, label, local_dir):
        Path(local_dir).mkdir(parents=True, exist_ok=True)
        out = Path(local_dir) / f"{label}.log"
        text = self.log_provider(label) if self.log_provider else self.log_text
        out.write_text(text)
        self.calls.append(("fetch_log", remote_job_id))
        return str(out)

    def cancel(self, remote_job_id):
        pass


def seed_running_job(conn, *, smiles="c1ccccc1", sig="sig-1", ssh_jobid="j-1"):
    return db.insert_job(
        conn,
        db.Job(
            id=None,
            signature=sig,
            smiles=smiles,
            method="b3lyp",
            basis="6-31g*",
            charge=0,
            mult=1,
            job_kind="properties",
            status="running",
            submitted_at="2026-05-21T00:00:00",
            started_at="2026-05-21T00:00:01",
            ssh_jobid=ssh_jobid,
        ),
    )
