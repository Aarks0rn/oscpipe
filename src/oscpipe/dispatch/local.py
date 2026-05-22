"""LocalBackend — runs g16 from $PATH in-place.

For tests + dev only. Production workflows go through SshBackend.
Submit is synchronous: it invokes g16 inline and writes the .log next to
the .com. jobid == label.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from typing import Literal

from ..dft.gaussian import is_log_complete

JobStatus = Literal["pending", "running", "complete", "error", "unknown"]


class LocalBackend:
    def __init__(self, work_dir: str, exe: str = "g16"):
        self.work_dir = Path(work_dir)
        self.exe = exe
        self.work_dir.mkdir(parents=True, exist_ok=True)

    def submit(self, com_path: str, label: str) -> str:
        com_src = Path(com_path)
        com_dest = self.work_dir / f"{label}.com"
        if com_src.resolve() != com_dest.resolve():
            shutil.copy(com_src, com_dest)
        # g16 writes <label>.log next to the .com itself; we only invoke it.
        subprocess.run(
            [self.exe, str(com_dest)],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            cwd=self.work_dir,
            check=False,
        )
        return label

    def poll(self, remote_job_id: str) -> JobStatus:
        log = self.work_dir / f"{remote_job_id}.log"
        if not log.exists():
            return "unknown"
        if is_log_complete(str(log)):
            return "complete"
        # If g16 returned but log lacks the marker → error.
        tail = log.read_bytes()[-4096:].decode("utf-8", errors="ignore")
        if "Error termination" in tail:
            return "error"
        return "running"

    def fetch_log(self, remote_job_id: str, label: str, local_dir: str) -> str:
        # Log is already local. If local_dir differs, copy it there.
        src = self.work_dir / f"{label}.log"
        os.makedirs(local_dir, exist_ok=True)
        dst = Path(local_dir) / f"{label}.log"
        if src.resolve() != dst.resolve():
            shutil.copy(src, dst)
        return str(dst)

    def cancel(self, remote_job_id: str) -> None:
        # Synchronous backend — nothing to cancel.
        pass

    def preflight(self) -> list[tuple[str, bool, str]]:
        """Health-check the local backend. Returns [(name, passed, message)]."""
        results = []

        exe_path = shutil.which(self.exe)
        results.append(("g16", bool(exe_path), exe_path or f"{self.exe!r} not in PATH"))

        probe = self.work_dir / ".preflight_probe"
        try:
            probe.write_text("x")
            probe.unlink()
            results.append(("work_dir", True, str(self.work_dir)))
        except OSError as exc:
            results.append(("work_dir", False, f"not writable: {exc}"))

        return results
