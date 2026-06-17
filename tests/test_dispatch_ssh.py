"""SshBackend connection + policy tests with paramiko mocked.

We inject a fake connect callable that yields a FakeSSH + FakeSFTP. The
paramiko module is never imported; tests run without network access. The
submit/poll/cancel/preflight launch lifecycle is covered in
test_dispatch_ssh_direct.py, which reuses the fakes defined here.
"""

from __future__ import annotations

import io
from pathlib import Path

import pytest

from oscpipe.dispatch.ssh import SshBackend
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


# ── fetch_log ──────────────────────────────────────────────────────────────


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


def test_fetch_log_downloads_to_local_dir(fake_settings, tmp_path):
    backend, _, sftp = _submitted_backend(fake_settings)
    out = tmp_path / "out"
    local = backend.fetch_log("12345", "x", str(out))
    assert local == str(out / "x.log")
    assert Path(local).read_text() == "downloaded log content\n"
    assert sftp.get_calls == [("/home/user/work/x.log", str(out / "x.log"))]


# ── reconnect on dropped connection ─────────────────────────────────────────


def test_run_reconnects_after_connection_reset(fake_settings):
    """A long λ_h run polls over one connection for hours; the workstation can
    reset it mid-run. _run must reconnect once and retry, not crash."""
    bad = FakeSSH()

    def boom(cmd):
        bad.cmds.append(cmd)
        raise ConnectionResetError(104, "Connection reset by peer")

    bad.exec_command = boom
    good = FakeSSH()
    good.default_stdout = b"ALIVE\n"
    conns = iter([(bad, FakeSFTP()), (good, FakeSFTP())])
    backend = SshBackend(fake_settings, _connect=lambda: next(conns))

    out, _ = backend._run("kill -0 1")
    assert "ALIVE" in out             # served by the healthy reconnect
    assert bad.cmds and good.cmds     # first conn tried, then retried on the second


def test_run_reconnects_after_ssh_session_not_active(fake_settings):
    """After hours of polling the paramiko transport can go inactive; exec_command
    then raises SSHException, not an OSError. _run must still reconnect+retry — this
    is the exact drop that crashed the diF sweep ('SSH session not active')."""
    from paramiko.ssh_exception import SSHException

    bad = FakeSSH()

    def boom(cmd):
        bad.cmds.append(cmd)
        raise SSHException("SSH session not active")

    bad.exec_command = boom
    good = FakeSSH()
    good.default_stdout = b"ALIVE\n"
    conns = iter([(bad, FakeSFTP()), (good, FakeSFTP())])
    backend = SshBackend(fake_settings, _connect=lambda: next(conns))

    out, _ = backend._run("kill -0 1")
    assert "ALIVE" in out             # served by the healthy reconnect
    assert bad.cmds and good.cmds     # first conn tried, then retried on the second


def test_run_reraises_if_reconnect_also_fails(fake_settings):
    """Retry is exactly once — a persistent failure propagates, no infinite loop."""

    def make_boom():
        s = FakeSSH()
        s.exec_command = lambda cmd: (_ for _ in ()).throw(ConnectionResetError(104, "reset"))
        return s, FakeSFTP()

    backend = SshBackend(fake_settings, _connect=lambda: make_boom())
    with pytest.raises(ConnectionResetError):
        backend._run("kill -0 1")
