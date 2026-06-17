"""Settings.load tests — env + kwargs precedence."""

from oscpipe.settings import Settings, load


def test_defaults():
    s = load()
    assert isinstance(s, Settings)
    assert s.gaussian_exe == "g16"


def test_env_override(monkeypatch):
    monkeypatch.setenv("OSC_REMOTE_HOST", "203.0.113.10")
    monkeypatch.setenv("OSC_GAUSSIAN_NPROC", "8")
    s = load()
    assert s.remote_host == "203.0.113.10"
    assert s.gaussian_nproc == 8


def test_kwargs_beat_env(monkeypatch):
    monkeypatch.setenv("OSC_REMOTE_HOST", "from-env")
    s = load(remote_host="from-kwargs")
    assert s.remote_host == "from-kwargs"


def test_max_lanes_default_is_one():
    # Dataclass default (not load(), which reads config.local.toml) — serial unless set.
    assert Settings().max_lanes == 1


def test_max_lanes_env_override(monkeypatch):
    monkeypatch.setenv("OSC_MAX_LANES", "2")
    s = load()
    assert s.max_lanes == 2
