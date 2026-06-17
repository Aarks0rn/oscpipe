"""SSH backend. paramiko + SFTP, strict host keys, key-only auth.

Security policy (locked):
    - paramiko.RejectPolicy + known_hosts for host keys
    - SSH private key only (settings.remote_key_file); password mode is disallowed
    - First-time setup: user runs `ssh user@host` once to pin the
      fingerprint, then `ssh-copy-id` to install the laptop's public key.

Implements the backend contract (submit / poll / fetch_log / cancel). g16 is launched
directly over SSH with `nohup ... & echo $!`; the shell PID is the job id, polled with
`kill -0`. The host has no batch scheduler (no qsub/qstat) — jobs run as soon as they
launch, so pace submissions yourself.
"""

from __future__ import annotations

import os
import re

from paramiko.ssh_exception import SSHException

from ..settings import Settings
from . import JobStatus


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
        # λ_h workflows poll over this one connection for ~a day; a keepalive
        # stops the workstation/NAT from dropping it during idle poll gaps.
        transport = client.get_transport()
        if transport is not None:
            transport.set_keepalive(30)
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

    def _reset(self) -> None:
        """Drop the (likely dead) connection so the next call reconnects."""
        for handle in (self._sftp, self._ssh):
            try:
                if handle is not None:
                    handle.close()
            except Exception:
                pass
        self._ssh = self._sftp = None

    # ── shell helper ───────────────────────────────────────────────────────

    def _exec(self, cmd: str) -> tuple[str, str]:
        self._ensure_open()
        _, stdout, stderr = self._ssh.exec_command(cmd)
        return stdout.read().decode(), stderr.read().decode()

    def _run(self, cmd: str, *, retry: bool = True) -> tuple[str, str]:
        try:
            return self._exec(cmd)
        except (OSError, EOFError, SSHException):
            # The workstation dropped an idle connection mid-workflow — seen as
            # ConnectionResetError (OSError) during a poll, or once the paramiko
            # transport has gone inactive, SSHException("SSH session not active")
            # (crashed the diF sweep, job 178). Reconnect once and retry so a
            # transient drop doesn't crash the whole multi-hour run.
            # retry=False for non-idempotent commands (launch): a reset between the
            # command running server-side and reading its reply would otherwise
            # re-submit and orphan the first job.
            if not retry:
                raise
            self._reset()
            return self._exec(cmd)

    # ── Backend protocol ───────────────────────────────────────────────────

    def submit(self, com_path: str, label: str) -> str:
        self._ensure_open()
        remote_dir = self.s.remote_work_dir.rstrip("/")
        _ensure_remote_dir(self._sftp, remote_dir)

        remote_com = f"{remote_dir}/{label}.com"
        remote_log = f"{remote_dir}/{label}.log"

        self._sftp.put(com_path, remote_com)

        # Clear any stale .log/.chk from a prior run of THIS label (the label hash
        # is route-independent, so a re-run reuses the filename). Without this the
        # poll() below could read the previous run's "Error termination".
        self._run(f"rm -f {remote_log} {remote_dir}/{label}.chk", retry=False)

        scrdir = f"GAUSS_SCRDIR={self.s.scratch_dir} " if self.s.scratch_dir else ""
        # nohup detaches g16 from the SSH channel so it survives once paramiko
        # closes the exec channel; `echo $!` prints the launched PID. retry=False:
        # a reconnect mid-launch would relaunch and orphan the first g16.
        launch = (
            f"cd {remote_dir} && "
            f"nohup {scrdir}{self.s.gaussian_exe} {label}.com {label}.log "
            f">/dev/null 2>&1 & echo $!"
        )
        out, err = self._run(launch, retry=False)
        pid = _parse_pid(out)
        if pid is None:
            raise RuntimeError(f"launch returned no PID. stdout={out!r} stderr={err!r}")
        self._jobs[pid] = remote_log
        return pid

    def poll(self, remote_job_id: str) -> JobStatus:
        out, _ = self._run(f"kill -0 {remote_job_id} 2>/dev/null && echo ALIVE || echo DEAD")
        if "ALIVE" in out:
            return "running"
        # Process is gone — decide from the log.
        log_path = self._jobs.get(remote_job_id)
        if log_path is None:
            return "unknown"
        status = self._poll_log(log_path)
        # _poll_log assumes a scheduler may not have started the job yet, so it
        # returns "pending"/"running" when the log is missing or unfinished. But
        # here the process has already exited, so an unfinished/absent log means
        # g16 died without completing → error.
        if status in ("pending", "running"):
            return "error"
        return status

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
        self._run(f"kill {remote_job_id} 2>/dev/null || true")

    def preflight(self) -> list[tuple[str, bool, str]]:
        """Health-check the host. Returns [(name, passed, message)]."""
        results = []

        # 1. g16 on PATH
        out, _ = self._run("which g16 2>/dev/null")
        g16_path = out.strip()
        results.append(("g16", bool(g16_path), g16_path or "not found in PATH"))

        # 2. work_dir writable (no scheduler on this host → no qstat check)
        remote_dir = self.s.remote_work_dir.rstrip("/")
        probe = f"{remote_dir}/.preflight_probe"
        out, err = self._run(f"mkdir -p {remote_dir} && touch {probe} && rm {probe} && echo ok")
        writable = out.strip() == "ok"
        results.append(
            ("work_dir", writable, remote_dir if writable else f"not writable: {err.strip()[:60]}")
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


def _parse_pid(stdout: str) -> str | None:
    """Return the last integer token in stdout (the `echo $!` PID), or None."""
    ids = re.findall(r"\d+", stdout)
    return ids[-1] if ids else None
