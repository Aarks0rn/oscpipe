"""SshBackend tests with paramiko mocked.

We inject a fake connect callable that yields a FakeSSH + FakeSFTP. The
paramiko module is never imported; tests run without network access.
"""

from __future__ import annotations

import io
from pathlib import Path

import pytest

from oscpipe.dispatch.ssh import SshBackend, _parse_qsub_jobid
from oscpipe.settings import Settings

# ── fakes ──────────────────────────────────────────────────────────────────


class _FakeFile:
    def __init__(self):
        self.buf = io.StringIO()

    def write(self, s):
        self.buf.write(s)

    def __enter__(self):
        return self

    def __exit__(self, *a):
        pass


class FakeSFTP:
    def __init__(self):
        self.put_calls: list[tuple[str, str]] = []
        self.get_calls: list[tuple[str, str]] = []
        self.opened: dict[str, _FakeFile] = {}
        self.existing_dirs: set[str] = {"/"}

    def stat(self, path):
        if path in self.existing_dirs:
            return object()
        raise FileNotFoundError(path)

    def mkdir(self, path):
        self.existing_dirs.add(path)

    def put(self, local, remote):
        self.put_calls.append((local, remote))

    def get(self, remote, local):
        Path(local).write_text("downloaded log content\n")
        self.get_calls.append((remote, local))

    def open(self, path, mode):
        f = _FakeFile()
        self.opened[path] = f
        return f

    def chmod(self, path, mode):
        pass

    def close(self):
        pass


class FakeStream:
    def __init__(self, data: bytes):
        self.data = data

    def read(self):
        return self.data


class FakeSSH:
    def __init__(self):
        self.cmds: list[str] = []
        # Default canned responses; tests can override.
        self.stdout_for: dict[str, bytes] = {}
        self.default_stdout = b""

    def exec_command(self, cmd):
        self.cmds.append(cmd)
        for needle, out in self.stdout_for.items():
            if needle in cmd:
                return None, FakeStream(out), FakeStream(b"")
        return None, FakeStream(self.default_stdout), FakeStream(b"")

    def close(self):
        pass


@pytest.fixture
def fake_settings(tmp_path):
    keyf = tmp_path / "id_rsa"
    keyf.write_text("-----BEGIN OPENSSH PRIVATE KEY-----\nfake\n")
    return Settings(
        remote_host="203.0.113.10",
        remote_user="user",
        remote_key_file=str(keyf),
        remote_work_dir="/home/user/work",
        gaussian_nproc=12,
    )


def _backend(settings, ssh: FakeSSH, sftp: FakeSFTP) -> SshBackend:
    return SshBackend(settings, _connect=lambda: (ssh, sftp))


# ── policy ─────────────────────────────────────────────────────────────────


def test_refuses_empty_key_file():
    s = Settings(remote_host="h", remote_user="u", remote_work_dir="/w")
    with pytest.raises(ValueError, match="key_file"):
        SshBackend(s)


def test_refuses_empty_host(tmp_path):
    keyf = tmp_path / "k"
    keyf.write_text("x")
    s = Settings(remote_key_file=str(keyf), remote_work_dir="/w")
    with pytest.raises(ValueError):
        SshBackend(s)


# ── submit ─────────────────────────────────────────────────────────────────


def test_submit_uploads_com_and_returns_jobid(fake_settings, tmp_path):
    com = tmp_path / "h2.com"
    com.write_text("%chk=h2.chk\n#p b3lyp/6-31g* opt\n\nh2\n\n0 1\nH 0 0 0\nH 0.74 0 0\n\n")

    ssh = FakeSSH()
    ssh.stdout_for["qsub"] = b'Your job 70867 ("h2") has been submitted\n'
    sftp = FakeSFTP()
    backend = _backend(fake_settings, ssh, sftp)

    jobid = backend.submit(str(com), "h2")
    assert jobid == "70867"
    # com was uploaded to the remote work dir
    assert any(local == str(com) and remote.endswith("/h2.com") for local, remote in sftp.put_calls)
    # qsub script was written
    assert any(path.endswith("/h2.qsub.sh") for path in sftp.opened)
    # qsub was run
    assert any("qsub h2.qsub.sh" in c for c in ssh.cmds)


def test_submit_raises_when_qsub_returns_no_jobid(fake_settings, tmp_path):
    com = tmp_path / "x.com"
    com.write_text("")
    ssh = FakeSSH()
    ssh.stdout_for["qsub"] = b"submit failed: queue is full\n"
    sftp = FakeSFTP()
    backend = _backend(fake_settings, ssh, sftp)
    with pytest.raises(RuntimeError, match="qsub"):
        backend.submit(str(com), "x")


# ── poll ───────────────────────────────────────────────────────────────────


def _submitted_backend(
    fake_settings, *, log_exists: bool = True
) -> tuple[SshBackend, FakeSSH, FakeSFTP]:
    ssh = FakeSSH()
    sftp = FakeSFTP()
    backend = _backend(fake_settings, ssh, sftp)
    # Bypass submit; seed the job→log map directly.
    log_path = "/home/user/work/x.log"
    backend._jobs["12345"] = log_path
    if log_exists:
        sftp.existing_dirs.add(log_path)
    return backend, ssh, sftp


def test_poll_running(fake_settings):
    backend, ssh, _ = _submitted_backend(fake_settings)
    ssh.stdout_for["qstat"] = b"job_number: 12345\njob_state: r running\n"
    assert backend.poll("12345") == "running"


def test_poll_complete(fake_settings):
    backend, ssh, _ = _submitted_backend(fake_settings)
    ssh.stdout_for["qstat"] = b"Following jobs do not exist: 12345\n"
    ssh.stdout_for["tail -c"] = b" Normal termination of Gaussian 16 at ...\n"
    assert backend.poll("12345") == "complete"


def test_poll_error(fake_settings):
    backend, ssh, _ = _submitted_backend(fake_settings)
    ssh.stdout_for["qstat"] = b"Following jobs do not exist: 12345\n"
    ssh.stdout_for["tail -c"] = b"Error termination via Lnk1e\n"
    assert backend.poll("12345") == "error"


def test_poll_eqw_returns_error(fake_settings):
    backend, ssh, _ = _submitted_backend(fake_settings)
    ssh.stdout_for["qstat"] = b"job 12345 state: Eqw\n"
    assert backend.poll("12345") == "error"


def test_poll_unknown_jobid(fake_settings):
    backend, ssh, _ = _submitted_backend(fake_settings)
    # No map entry for this jobid.
    backend._jobs.clear()
    ssh.stdout_for["qstat"] = b"Following jobs do not exist: 99\n"
    assert backend.poll("99") == "unknown"


def test_poll_pending_when_log_absent(fake_settings):
    """UGE accounting race: qstat says 'do not exist' before scheduler picks it up,
    log file not yet created — must be 'pending', not 'unknown' (would mark job error)."""
    backend, ssh, _ = _submitted_backend(fake_settings, log_exists=False)
    ssh.stdout_for["qstat"] = b"Following jobs do not exist: 12345\n"
    assert backend.poll("12345") == "pending"


# ── fetch_log ──────────────────────────────────────────────────────────────


def test_fetch_log_downloads_to_local_dir(fake_settings, tmp_path):
    backend, _, sftp = _submitted_backend(fake_settings)
    out = tmp_path / "out"
    local = backend.fetch_log("12345", "x", str(out))
    assert local == str(out / "x.log")
    assert Path(local).read_text() == "downloaded log content\n"
    assert sftp.get_calls == [("/home/user/work/x.log", str(out / "x.log"))]


# ── cancel ─────────────────────────────────────────────────────────────────


def test_cancel_runs_qdel(fake_settings):
    backend, ssh, _ = _submitted_backend(fake_settings)
    backend.cancel("12345")
    assert any("qdel 12345" in c for c in ssh.cmds)


# ── qsub parser unit ───────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "stdout,expected",
    [
        ('Your job 70867 ("h2") has been submitted\n', "70867"),
        ("Your job-array 1234.1-10:1 has been submitted\n", "1234"),
        ("submit failed: queue full\n", None),
    ],
)
def test_parse_qsub_jobid(stdout, expected):
    assert _parse_qsub_jobid(stdout) == expected
