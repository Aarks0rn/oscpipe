"""SshBackend launch-lifecycle tests — direct g16 launch over SSH.

Reuses the paramiko fakes from test_dispatch_ssh. The backend launches
g16 with `nohup ... & echo $!`, tracks the shell PID, and polls with `kill -0`.
"""

from __future__ import annotations

import pytest
from test_dispatch_ssh import FakeSFTP, FakeSSH

from oscpipe.dispatch.ssh import SshBackend, _parse_pid
from oscpipe.settings import Settings


@pytest.fixture
def fake_settings(tmp_path):
    keyf = tmp_path / "id_rsa"
    keyf.write_text("-----BEGIN OPENSSH PRIVATE KEY-----\nfake\n")
    return Settings(
        remote_host="10.93.48.2",
        remote_user="MiniLab",
        remote_key_file=str(keyf),
        remote_work_dir="/home/MiniLab/Lab-Data/SA26/oscpipe-work",
        gaussian_nproc=12,
    )


def _backend(settings, ssh: FakeSSH, sftp: FakeSFTP) -> SshBackend:
    return SshBackend(settings, _connect=lambda: (ssh, sftp))


def _submit(backend, ssh, sftp, com, label="job_x"):
    return backend.submit(str(com), label)


# ── PID parsing ─────────────────────────────────────────────────────────────


def test_parse_pid_takes_last_integer():
    assert _parse_pid("12345\n") == "12345"
    assert _parse_pid("some noise\n67890\n") == "67890"
    assert _parse_pid("no pid here") is None


# ── submit ──────────────────────────────────────────────────────────────────


def test_submit_launches_nohup_g16_and_returns_pid(fake_settings, tmp_path):
    ssh, sftp = FakeSSH(), FakeSFTP()
    ssh.stdout_for = {"echo $!": b"4242\n"}
    com = tmp_path / "job_x.com"
    com.write_text("# header\n")
    backend = _backend(fake_settings, ssh, sftp)

    jobid = _submit(backend, ssh, sftp, com)

    assert jobid == "4242"
    # .com uploaded
    assert any(remote.endswith("job_x.com") for _, remote in sftp.put_calls)
    # a launch command ran g16 in the background
    launch = [c for c in ssh.cmds if "nohup" in c]
    assert launch, ssh.cmds
    assert "job_x.com" in launch[0] and "job_x.log" in launch[0]
    assert "echo $!" in launch[0]
    # jobid maps to the remote log so poll/fetch can find it
    assert backend._jobs["4242"].endswith("job_x.log")
    # NO qsub script written
    assert not any("qsub" in c for c in ssh.cmds)


def test_submit_clears_stale_log_before_launch(fake_settings, tmp_path):
    ssh, sftp = FakeSSH(), FakeSFTP()
    ssh.stdout_for = {"echo $!": b"9\n"}
    com = tmp_path / "job_x.com"
    com.write_text("# header\n")
    backend = _backend(fake_settings, ssh, sftp)

    _submit(backend, ssh, sftp, com)

    rm = [c for c in ssh.cmds if c.startswith("rm -f") or "rm -f" in c]
    assert any("job_x.log" in c and "job_x.chk" in c for c in rm), ssh.cmds


# ── poll ────────────────────────────────────────────────────────────────────


def test_poll_running_when_process_alive(fake_settings):
    ssh, sftp = FakeSSH(), FakeSFTP()
    ssh.stdout_for = {"kill -0": b"ALIVE\n"}
    backend = _backend(fake_settings, ssh, sftp)
    backend._jobs["4242"] = "/home/MiniLab/Lab-Data/SA26/oscpipe-work/job_x.log"

    assert backend.poll("4242") == "running"


def test_poll_complete_when_dead_and_log_normal(fake_settings):
    ssh, sftp = FakeSSH(), FakeSFTP()
    log = "/home/MiniLab/Lab-Data/SA26/oscpipe-work/job_x.log"
    sftp.existing_dirs.add(log)  # log exists
    ssh.stdout_for = {
        "kill -0": b"DEAD\n",
        "tail -c": b"   Normal termination of Gaussian 16\n",
    }
    backend = _backend(fake_settings, ssh, sftp)
    backend._jobs["4242"] = log

    assert backend.poll("4242") == "complete"


def test_poll_error_when_dead_and_log_has_error(fake_settings):
    ssh, sftp = FakeSSH(), FakeSFTP()
    log = "/home/MiniLab/Lab-Data/SA26/oscpipe-work/job_x.log"
    sftp.existing_dirs.add(log)
    ssh.stdout_for = {
        "kill -0": b"DEAD\n",
        "tail -c": b" Error termination via Lnk1e\n",
    }
    backend = _backend(fake_settings, ssh, sftp)
    backend._jobs["4242"] = log

    assert backend.poll("4242") == "error"


def test_poll_error_when_dead_and_log_missing(fake_settings):
    """Process gone but g16 never wrote a finished log → it crashed → error."""
    ssh, sftp = FakeSSH(), FakeSFTP()
    log = "/home/MiniLab/Lab-Data/SA26/oscpipe-work/job_x.log"
    # log NOT in existing_dirs → sftp.stat raises → _poll_log would say 'pending'
    ssh.stdout_for = {"kill -0": b"DEAD\n"}
    backend = _backend(fake_settings, ssh, sftp)
    backend._jobs["4242"] = log

    assert backend.poll("4242") == "error"


def test_poll_unknown_when_pid_untracked(fake_settings):
    ssh, sftp = FakeSSH(), FakeSFTP()
    ssh.stdout_for = {"kill -0": b"DEAD\n"}
    backend = _backend(fake_settings, ssh, sftp)
    assert backend.poll("999") == "unknown"


# ── cancel ──────────────────────────────────────────────────────────────────


def test_cancel_kills_pid(fake_settings):
    ssh, sftp = FakeSSH(), FakeSFTP()
    backend = _backend(fake_settings, ssh, sftp)
    backend.cancel("4242")
    assert any("kill 4242" in c for c in ssh.cmds), ssh.cmds


# ── preflight ───────────────────────────────────────────────────────────────


def test_preflight_checks_g16_and_workdir_no_qstat(fake_settings):
    ssh, sftp = FakeSSH(), FakeSFTP()
    ssh.stdout_for = {
        "which g16": b"/usr/local/g16/g16\n",
        "echo ok": b"ok\n",
    }
    backend = _backend(fake_settings, ssh, sftp)

    checks = backend.preflight()
    names = {name for name, _, _ in checks}

    assert "g16" in names
    assert "work_dir" in names or "scratch" in names
    assert "qstat" not in names  # no scheduler on a direct host
    assert all(passed for _, passed, _ in checks), checks
