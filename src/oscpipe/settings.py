"""Immutable settings loader.

Single source of truth for compute and IO config. Loaded once at entry-point;
pass the Settings instance through; never read os.environ ad-hoc.

Precedence (highest first):
    1. Constructor kwargs (tests)
    2. Environment variables (OSC_*)
    3. config.toml next to the user's working dir
    4. Hardcoded defaults below.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Settings:
    remote_host: str = ""  # e.g. 203.0.113.10
    remote_user: str = ""  # e.g. user
    remote_port: int = 22
    remote_key_file: str = ""  # absolute path to SSH private key (no password mode)
    remote_work_dir: str = ""  # workstation directory for .com/.log
    known_hosts_path: str = ""  # absolute path; strict policy when set

    gaussian_exe: str = "g16"
    gaussian_nproc: int = 12
    gaussian_mem: str = "28GB"
    scratch_dir: str = ""

    db_path: str = "results.db"
    project: str = ""

    # Dispatch
    backend: str = "ssh"  # ssh | local
    local_work_dir: str = ""  # used when backend=local
    poll_interval_seconds: int = 5
    # Max g16 jobs run concurrently within one workflow layer. 1 = serial (every
    # existing single-call path is unchanged). On a direct host with no queue the
    # operator sets this; invariant: max_lanes * gaussian_nproc <= host cores and
    # max_lanes * gaussian_mem <~ host RAM - headroom.
    max_lanes: int = 1


_ENV_PREFIX = "OSC_"


def _load_toml() -> dict:
    """Return dict from config.local.toml or config.toml in CWD, or {} if neither exists."""
    import tomllib
    from pathlib import Path

    for name in ("config.local.toml", "config.toml"):
        p = Path(name)
        if p.exists():
            with p.open("rb") as f:
                return tomllib.load(f)
    return {}


def load(**overrides) -> Settings:
    """Load Settings from environment + kwargs.

    Precedence: kwargs > OSC_<FIELD> env vars > config.local.toml / config.toml > defaults.
    Unknown kwargs raise TypeError (via Settings constructor). Env values and TOML
    values are coerced to int for numeric fields; everything else is read as a string.
    """
    import os
    from dataclasses import fields

    toml_data = _load_toml()
    kwargs: dict = {}
    for f in fields(Settings):
        # env wins over toml; both coerce to int for int fields, else str.
        if (env_key := _ENV_PREFIX + f.name.upper()) in os.environ:
            raw = os.environ[env_key]
        elif f.name in toml_data:
            raw = toml_data[f.name]
        else:
            continue
        is_int = f.type is int or f.type == "int"
        kwargs[f.name] = int(raw) if is_int else str(raw)
    kwargs.update(overrides)
    return Settings(**kwargs)
