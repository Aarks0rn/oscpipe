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
    remote_pe: str = "OpenMP"  # UGE parallel environment name (e.g. OpenMP, mpi)

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


_ENV_PREFIX = "OSC_"


def load(**overrides) -> Settings:
    """Load Settings from environment + kwargs.

    Precedence: kwargs > OSC_<FIELD> env vars > dataclass defaults.
    Unknown kwargs raise TypeError (via Settings constructor). Env values are
    coerced to int for numeric fields; everything else is read as a string.
    """
    import os
    from dataclasses import fields

    kwargs: dict = {}
    for f in fields(Settings):
        env_key = _ENV_PREFIX + f.name.upper()
        if env_key in os.environ:
            raw = os.environ[env_key]
            kwargs[f.name] = int(raw) if f.type is int or f.type == "int" else raw
    kwargs.update(overrides)
    return Settings(**kwargs)
