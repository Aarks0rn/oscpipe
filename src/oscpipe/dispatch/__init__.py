"""Dispatch backends. JobStatus is the shared job-state alias every backend returns."""

from __future__ import annotations

from typing import Literal

JobStatus = Literal["pending", "running", "complete", "error", "unknown"]
