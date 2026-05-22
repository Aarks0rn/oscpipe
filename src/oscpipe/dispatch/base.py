"""Backend protocol — fire-and-forget submit + status poll + fetch.

Every backend (only `ssh` for first cut) implements this protocol so the
CLI and Streamlit dashboard never branch on backend.
"""

from __future__ import annotations

from typing import Literal, Protocol

JobStatus = Literal["pending", "running", "complete", "error", "unknown"]


class Backend(Protocol):
    def submit(self, com_path: str, label: str) -> str:
        """Upload .com to compute host, queue it, return remote job id."""
        ...

    def poll(self, remote_job_id: str) -> JobStatus:
        """Return current job state. Cheap call; safe to invoke per refresh."""
        ...

    def fetch_log(self, remote_job_id: str, label: str, local_dir: str) -> str:
        """Download the completed .log to local_dir; return local path."""
        ...

    def cancel(self, remote_job_id: str) -> None: ...
