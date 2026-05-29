"""Job runner — owns the lifecycle of a single Gaussian job.

A :class:`JobRunner` bundles the ``(conn, settings, backend)`` triple plus the
stdout sink and runs one ``db.Job`` through submit → poll/wait → fetch → parse →
persist.  The remote-label naming convention and the per-``job_kind`` result
persisting dispatch live here, not scattered across the CLI.

Layer: sits above ``store`` / ``dft`` / ``dispatch``.  ``cli`` and ``workflows``
are its callers; nothing here imports upward.
"""

from __future__ import annotations

import datetime
import json
import os
import tempfile
import time
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import IO, NamedTuple

from .dft import gaussian
from .settings import Settings
from .store import db


def _now() -> str:
    return datetime.datetime.now().isoformat(timespec="seconds")


# ── label + path conventions ─────────────────────────────────────────────────


def label(smiles_canonical: str, sig: str) -> str:
    """Remote .com/.log/.chk basename for a job: <slug>_<signature>."""
    slug = "".join(c if c.isalnum() else "_" for c in smiles_canonical)[:24]
    return f"{slug}_{sig}"


def label_for_row(row) -> str:
    """Rebuild a job's remote label from its DB row.

    Lambda sub-jobs store a plain identifier in ``notes`` ("neutral_opt") used
    as a label suffix.  Metadata notes ("nstates=10") contain "=" and must NOT
    be appended.
    """
    base = label(row["smiles"], row["signature"])
    if row["notes"] and "=" not in row["notes"]:
        return f"{base}_{row['notes']}"
    return base


def db_dir(settings: Settings) -> Path:
    return Path(settings.db_path).resolve().parent


def log_dir(settings: Settings) -> Path:
    d = db_dir(settings) / "logs"
    d.mkdir(parents=True, exist_ok=True)
    return d


def remote_log_path(settings: Settings, label_: str) -> str:
    base = (settings.remote_work_dir or "").rstrip("/")
    return f"{base}/{label_}.log" if base else f"{label_}.log"


class Resolved(NamedTuple):
    """Outcome of resolving one job.

    Unpacks as ``(job_id, status, log_path, result, cached)``.  ``result`` is a
    :class:`gaussian.PropertiesResult` for properties / sp_neutral / sp_cation
    jobs and ``None`` for every other kind; ``cached`` is True when the result
    came from the signature cache without re-running Gaussian.
    """

    job_id: int | None
    status: str
    log_path: str | None
    result: object | None
    cached: bool


@dataclass
class JobRunner:
    """Runs single Gaussian jobs through their full lifecycle.

    Bundles the connection, settings, backend and stdout so callers stop
    threading them through every helper.
    """

    conn: object
    settings: Settings
    backend: object
    stdout: IO[str]

    # ── submit ───────────────────────────────────────────────────────────────

    def _submit(self, job: db.Job, com_text: str, label_: str) -> tuple[int, str]:
        """Insert pending row, hand .com to backend, mark running.

        Returns (job_id, ssh_jobid).
        """
        job_id = db.insert_job(self.conn, job)
        with tempfile.TemporaryDirectory() as tmp:
            com_path = Path(tmp) / f"{label_}.com"
            com_path.write_text(com_text)
            ssh_jobid = self.backend.submit(str(com_path), label_)
        db.update_job_status(self.conn, job_id, "running", started_at=_now(), ssh_jobid=ssh_jobid)
        print(f"submitted: job={job_id} jobid={ssh_jobid} label={label_}", file=self.stdout)
        return job_id, ssh_jobid

    def resolve(
        self,
        job: db.Job,
        label_: str,
        build_com: Callable[[], str],
        *,
        wait: bool,
        need_log: bool = False,
    ) -> Resolved:
        """Resolve one job to its result: cache-hit, or submit → poll/wait → persist.

        On a signature cache hit the result is returned without touching the
        backend and ``build_com`` is never called (so the embed/write is skipped).
        On a miss the job is submitted, then polled once (``wait=False``,
        fire-and-forget) or blocked on until terminal (``wait=True``); a completed
        job is fetched, parsed and persisted by ``job.job_kind``.

        ``need_log=True`` invalidates a cache hit whose local .log is gone — the
        caller needs the geometry, so the job must be recomputed.
        """
        hit = db.find_complete_by_signature(self.conn, job.signature)
        if hit is not None:
            cached = self._load_cached(hit, need_log=need_log)
            if cached is not None:
                return cached

        job_id, ssh_jobid = self._submit(job, build_com(), label_)
        status = self._wait(ssh_jobid) if wait else self.backend.poll(ssh_jobid)

        if status == "detached":
            print(
                f"detached: job={job_id} still running — run `oscpipe fetch` to pull results",
                file=self.stdout,
            )
            return Resolved(job_id, status, None, None, False)
        if status != "complete":
            if status == "error":
                db.update_job_status(
                    self.conn,
                    job_id,
                    "error",
                    completed_at=_now(),
                    error_msg="backend reported error",
                )
                return Resolved(job_id, "error", None, None, False)
            return Resolved(job_id, status, None, None, False)  # async: pending/running

        local_log = self.backend.fetch_log(ssh_jobid, label_, str(log_dir(self.settings)))
        final, result = self._persist_by_kind(job_id, job.job_kind, local_log)
        return Resolved(job_id, final, local_log, result, False)

    def _load_cached(self, hit, *, need_log: bool) -> Resolved | None:
        """Rebuild a :class:`Resolved` from a cache-hit row, or None if unusable.

        Unusable = the stored result is incomplete (e.g. a pre-existing job that
        reached 'complete' without a results row) or ``need_log`` is set but the
        local .log no longer exists.
        """
        log_path = hit["log_path"]
        if need_log and (not log_path or not os.path.exists(log_path)):
            return None
        kind = hit["job_kind"]
        if kind == "tddft":
            spec = self.conn.execute(
                "SELECT spectra_json FROM results WHERE job_id = ?", (hit["id"],)
            ).fetchone()
            if not spec or not spec["spectra_json"]:
                return None
            return Resolved(hit["id"], "complete", log_path, None, True)
        if kind == "sp_dimer":
            return Resolved(hit["id"], "complete", log_path, None, True)
        # properties / sp_neutral / sp_cation — require a real results row. A
        # pre-refactor λ_h sub-job reached 'complete' with no row (all fields
        # NULL via the LEFT JOIN); recompute rather than feed NULLs downstream.
        # energy_ev may legitimately be absent on an externally-seeded row, so it
        # is not part of the guard (pipeline rows always carry it).
        if any(hit[c] is None for c in ("homo_ev", "lumo_ev", "gap_ev")):
            return None
        props = gaussian.PropertiesResult(
            homo_ev=hit["homo_ev"],
            lumo_ev=hit["lumo_ev"],
            gap_ev=hit["gap_ev"],
            dipole_debye=hit["dipole_debye"],
            energy_ev=hit["energy_ev"],
        )
        return Resolved(hit["id"], "complete", log_path, props, True)

    def _wait(self, ssh_jobid: str) -> str:
        """Block until the backend reports a terminal state.

        Ctrl-C detaches cleanly — the job keeps running on the workstation.
        """
        try:
            while True:
                s = self.backend.poll(ssh_jobid)
                if s in ("complete", "error", "unknown"):
                    return s
                time.sleep(self.settings.poll_interval_seconds)
        except KeyboardInterrupt:
            return "detached"

    # ── persist by job_kind (single home) ──────────────────────────────────────

    def _persist_by_kind(
        self, job_id: int, job_kind: str, local_log: str
    ) -> tuple[str, object | None]:
        """Parse a completed log by ``job_kind`` and persist it. Returns (status, result).

        The one place that maps a ``job_kind`` to how its result is parsed and
        stored — total over every kind the schema defines.  ``result`` is the
        PropertiesResult for properties / sp_neutral / sp_cation, else None.
        A parse failure marks the job 'error'; an unknown kind (e.g. ``freq``,
        never produced yet) raises rather than persisting silently.
        """
        persisters = {
            "properties": self._persist_properties,
            "sp_neutral": self._persist_properties,
            "sp_cation": self._persist_properties,
            "tddft": self._persist_tddft,
            "sp_dimer": self._persist_dimer,
        }
        persist = persisters.get(job_kind)
        if persist is None:
            raise ValueError(f"no persister for job_kind={job_kind!r}")
        try:
            return persist(job_id, local_log)
        except ValueError as exc:
            db.update_job_status(
                self.conn,
                job_id,
                "error",
                completed_at=_now(),
                log_path=local_log,
                error_msg=str(exc),
            )
            print(f"job {job_id} complete but parse failed: {exc}", file=self.stdout)
            return "error", None

    def _persist_properties(self, job_id: int, local_log: str) -> tuple[str, object]:
        props = gaussian.parse_properties(local_log)
        db.insert_result(
            self.conn,
            db.Result(
                job_id=job_id,
                homo_ev=props.homo_ev,
                lumo_ev=props.lumo_ev,
                gap_ev=props.gap_ev,
                dipole_debye=props.dipole_debye,
                energy_ev=props.energy_ev,
            ),
        )
        db.update_job_status(self.conn, job_id, "complete", completed_at=_now(), log_path=local_log)
        print(
            f"complete: job={job_id} HOMO={props.homo_ev:.3f} eV "
            f"LUMO={props.lumo_ev:.3f} eV gap={props.gap_ev:.3f} eV",
            file=self.stdout,
        )
        return "complete", props

    def _persist_tddft(self, job_id: int, local_log: str) -> tuple[str, None]:
        states = gaussian.parse_excited_states(local_log)
        if not states:
            raise ValueError(f"{local_log}: no excited states parsed")
        spectra = [
            {
                "n": s.n,
                "energy_ev": s.energy_ev,
                "wavelength_nm": s.wavelength_nm,
                "f": s.oscillator_strength,
            }
            for s in states
        ]
        db.insert_result(
            self.conn,
            db.Result(job_id=job_id, spectra_json=json.dumps(spectra)),
        )
        db.update_job_status(self.conn, job_id, "complete", completed_at=_now(), log_path=local_log)
        bright = max(states, key=lambda s: s.oscillator_strength)
        print(
            f"complete: job={job_id} tddft n={len(states)} "
            f"brightest λ={bright.wavelength_nm:.1f} nm f={bright.oscillator_strength:.3f}",
            file=self.stdout,
        )
        return "complete", None

    def _persist_dimer(self, job_id: int, local_log: str) -> tuple[str, None]:
        """Dimer SP keeps no results row — its J_hole is derived by the λ_h
        workflow (analysis.indo) and stored in the workflow summary, not per-job.
        Only the job status/log_path are recorded here."""
        db.update_job_status(self.conn, job_id, "complete", completed_at=_now(), log_path=local_log)
        return "complete", None

    # ── finish: resume an already-complete row (fetch / reconcile path) ─────────

    def finish_completed(self, row) -> str:
        """Fetch the log for a row already known complete, then persist by kind."""
        label_ = label_for_row(row)
        local_log = self.backend.fetch_log(row["ssh_jobid"], label_, str(log_dir(self.settings)))
        status, _ = self._persist_by_kind(row["id"], row["job_kind"], local_log)
        return status

    # ── rehydrate ──────────────────────────────────────────────────────────────

    def rehydrate(self, rows) -> None:
        """For SshBackend, repopulate jobid → remote_log_path from DB rows."""
        if not hasattr(self.backend, "_jobs"):
            return
        for r in rows:
            if not r["ssh_jobid"]:
                continue
            self.backend._jobs[r["ssh_jobid"]] = remote_log_path(self.settings, label_for_row(r))
