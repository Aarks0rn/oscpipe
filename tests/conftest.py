"""Pytest config — shared fixtures and marker registration."""

import os
from pathlib import Path

import pytest

HERE = Path(__file__).parent


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "real: requires live workstation (user@203.0.113.10) with g16 — run with -m real",
    )


def pytest_collection_modifyitems(config, items):
    """Auto-skip @pytest.mark.real tests unless -m real was explicitly requested."""
    if "real" in (config.option.markexpr or ""):
        return
    skip = pytest.mark.skip(reason="workstation test — run with -m real")
    for item in items:
        if item.get_closest_marker("real"):
            item.add_marker(skip)


@pytest.fixture
def fixtures_dir():
    return HERE / "fixtures"


@pytest.fixture
def real_settings(tmp_path):
    """Settings pointing at the live workstation. Skips if no SSH key found."""
    from oscpipe.settings import Settings

    key_file = os.environ.get("OSC_REMOTE_KEY_FILE", "")
    if not key_file:
        for candidate in ("~/.ssh/id_ed25519", "~/.ssh/id_rsa"):
            if Path(os.path.expanduser(candidate)).exists():
                key_file = candidate
                break
    if not key_file:
        pytest.skip("no SSH key found — set OSC_REMOTE_KEY_FILE or run ssh-keygen")

    return Settings(
        backend="ssh",
        remote_host="203.0.113.10",
        remote_user="user",
        remote_port=22,
        remote_key_file=key_file,
        remote_work_dir="/home/user/Lab-Data/SA26",
        known_hosts_path="~/.ssh/known_hosts",
        gaussian_exe="g16",
        gaussian_nproc=1,
        gaussian_mem="1GB",
        db_path=str(tmp_path / "results.db"),
        poll_interval_seconds=5,
    )
