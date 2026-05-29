"""SSH backend. paramiko + SFTP, strict host keys, key-only auth.

Security policy (locked):
    - paramiko.RejectPolicy + known_hosts for host keys
    - SSH private key only (settings.remote_key_file); password mode is disallowed
    - First-time setup: user runs `ssh user@203.0.113.10` once to pin the
      fingerprint, then `ssh-copy-id` to install the laptop's public key.

Implements `oscpipe.dispatch.base.Backend`. UGE (qsub) on the workstation.
"""

from __future__ import annotations

import os
import re
from typing import Literal

from ..settings import Settings

JobStatus = Literal["pending", "running", "complete", "error", "unknown"]


class SshBackend:
    def __init__(self, settings: Settings, *, _connect=None):
        if not settings.remote_key_file:
            raise ValueError(
                "remote_key_file is required; password auth is disallowed. "
                "See docs/ARCHITECTURE.md for the SSH setup steps."
            )
        if not settings.remote_host or not settings.remote_user:
            raise ValueError("remote_host and remote_user must be set")
        if not settings.remote_work_dir:
            raise ValueError("remote_work_dir must be set")
        self.s = settings
        self._connect = _connect or self._paramiko_connect
        self._ssh = None
        self._sftp = None
        # {remote_job_id: remote_log_path} so poll/fetch can locate the log.
        self._jobs: dict[str, str] = {}

    # ── connection ─────────────────────────────────────────────────────────

    def _paramiko_connect(self):
        import paramiko

        client = paramiko.SSHClient()
        if self.s.known_hosts_path:
            client.load_host_keys(os.path.expanduser(self.s.known_hosts_path))
        client.set_missing_host_key_policy(paramiko.RejectPolicy())
        client.connect(
            hostname=self.s.remote_host,
            port=self.s.remote_port,
            username=self.s.remote_user,
            key_filename=os.path.expanduser(self.s.remote_key_file),
            timeout=10,
            look_for_keys=False,
            allow_agent=False,
        )
        return client, client.open_sftp()

    def _ensure_open(self):
        if self._ssh is None:
            self._ssh, self._sftp = self._connect()

    def close(self) -> None:
        if self._sftp is not None:
            self._sftp.close()
        if self._ssh is not None:
            self._ssh.close()
        self._ssh = self._sftp = None

    # ── shell helper ───────────────────────────────────────────────────────

    def _run(self, cmd: str) -> tuple[str, str]:
        self._ensure_open()
        _, stdout, stderr = self._ssh.exec_command(cmd)
        return stdout.read().decode(), stderr.read().decode()

    # ── Backend protocol ───────────────────────────────────────────────────

    def submit(self, com_path: str, label: str) -> str:
        self._ensure_open()
        remote_dir = self.s.remote_work_dir.rstrip("/")
        _ensure_remote_dir(self._sftp, remote_dir)

        remote_com = f"{remote_dir}/{label}.com"
        remote_log = f"{remote_dir}/{label}.log"
        remote_script = f"{remote_dir}/{label}.qsub.sh"

        self._sftp.put(com_path, remote_com)

        script = _qsub_script(
            remote_dir=remote_dir,
            label=label,
            nproc=self.s.gaussian_nproc,
            pe=self.s.remote_pe,
            exe=self.s.gaussian_exe,
            scratch_dir=self.s.scratch_dir,
        )
        # paramiko exec for the heredoc-free write: use a sftp file handle.
        with self._sftp.open(remote_script, "w") as f:
            f.write(script)
        self._sftp.chmod(remote_script, 0o755)

        out, err = self._run(f"cd {remote_dir} && qsub {label}.qsub.sh")
        jobid = _parse_qsub_jobid(out)
        if jobid is None:
            raise RuntimeError(f"qsub returned no job id. stdout={out!r} stderr={err!r}")
        self._jobs[jobid] = remote_log
        return jobid

    def poll(self, remote_job_id: str) -> JobStatus:
        out, _ = self._run(f"qstat -j {remote_job_id} 2>&1")
        low = out.lower()
        if "eqw" in low or "error" in low and "do not exist" not in low:
            # qstat reports a queue-side error.
            return "error"
        if "do not exist" in low or "not found" in low:
            # Job left the queue. Inspect the log to decide complete vs error.
            log_path = self._jobs.get(remote_job_id)
            if log_path is None:
                return "unknown"
            return self._poll_log(log_path)
        if " r " in low or "running" in low:
            return "running"
        return "pending"

    def _poll_log(self, remote_log_path: str) -> JobStatus:
        # If the log doesn't exist yet, UGE hasn't started the job — treat as pending,
        # not unknown. Avoids a race right after qsub where qstat says "do not exist"
        # before the scheduler has populated job state.
        try:
            self._ensure_open()
            self._sftp.stat(remote_log_path)
        except IOError:
            return "pending"
        tail, _ = self._run(f"tail -c 4096 {remote_log_path} 2>/dev/null")
        if not tail:
            return "running"
        if "Normal termination of Gaussian" in tail:
            return "complete"
        if "Error termination" in tail:
            return "error"
        return "running"

    def fetch_log(self, remote_job_id: str, label: str, local_dir: str) -> str:
        self._ensure_open()
        remote_log = self._jobs.get(
            remote_job_id, f"{self.s.remote_work_dir.rstrip('/')}/{label}.log"
        )
        os.makedirs(local_dir, exist_ok=True)
        local_path = os.path.join(local_dir, f"{label}.log")
        self._sftp.get(remote_log, local_path)
        return local_path

    def cancel(self, remote_job_id: str) -> None:
        self._run(f"qdel {remote_job_id}")

    def preflight(self) -> list[tuple[str, bool, str]]:
        """Health-check the remote workstation. Returns [(name, passed, message)]."""
        results = []

        # 1. g16 on PATH
        out, _ = self._run("which g16 2>/dev/null")
        g16_path = out.strip()
        results.append(("g16", bool(g16_path), g16_path or "not found in PATH"))

        # 2. qstat exists (job 0 won't exist; UGE still prints "do not exist: 0")
        out, _ = self._run("qstat -j 0 2>&1 || true")
        low = out.lower()
        qstat_ok = "do not exist" in low or "unknown job" in low or "following" in low
        results.append(("qstat", qstat_ok, out.strip()[:80] if qstat_ok else "qstat not found"))

        # 3. scratch_dir writable
        remote_dir = self.s.remote_work_dir.rstrip("/")
        probe = f"{remote_dir}/.preflight_probe"
        out, err = self._run(f"touch {probe} && rm {probe} && echo ok")
        writable = out.strip() == "ok"
        results.append(
            ("scratch", writable, remote_dir if writable else f"not writable: {err.strip()[:60]}")
        )

        return results


# ── module helpers ─────────────────────────────────────────────────────────


def _ensure_remote_dir(sftp, remote_dir: str) -> None:
    parts = remote_dir.strip("/").split("/")
    path = ""
    for part in parts:
        path += "/" + part
        try:
            sftp.stat(path)
        except FileNotFoundError:
            sftp.mkdir(path)


def _qsub_script(
    *, remote_dir: str, label: str, nproc: int, pe: str, exe: str, scratch_dir: str = ""
) -> str:
    pe_line = f"#$ -pe {pe} {nproc}\n" if nproc > 1 else ""
    job_name = label[:64]  # UGE job name (max ~64 chars)
    scrdir_line = f"export GAUSS_SCRDIR={scratch_dir}\n" if scratch_dir else ""
    return (
        "#!/bin/sh\n"
        "#$ -S /bin/sh\n"
        "#$ -V\n"
        f"#$ -N {job_name}\n" + pe_line + "#$ -q all.q\n"
        "#$ -cwd\n"
        "#$ -j y\n"
        f"cd {remote_dir}\n" + scrdir_line + f"{exe} {label}.com {label}.log\n"
    )


def _parse_qsub_jobid(qsub_stdout: str) -> str | None:
    m = re.search(r"\bjob[- ]?(?:array )?(\d+)", qsub_stdout, re.IGNORECASE)
    return m.group(1) if m else None
