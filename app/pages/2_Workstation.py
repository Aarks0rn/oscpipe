"""Workstation page — live job dashboard, fetch, reconcile, auto-refresh."""

from __future__ import annotations

import io
import time
import types

import pandas as pd
import streamlit as st

from oscpipe.cli.main import _make_backend, run_fetch, run_reconcile
from oscpipe.settings import load
from oscpipe.store import db

st.title("Workstation")

settings = load()
host_label = f"{settings.remote_user}@{settings.remote_host}" if settings.remote_host else "local"
st.caption(f"Backend: {settings.backend}  |  {host_label}")


def _open_db():
    return db.open(settings.db_path)


def _try_backend():
    try:
        return _make_backend(settings), None
    except ValueError as e:
        return None, str(e)


# ── action bar ────────────────────────────────────────────────────────────────
btn_fetch, btn_reconcile, _, toggle_col = st.columns([1, 1, 3, 2])

do_fetch = btn_fetch.button("↺ Fetch")
do_reconcile = btn_reconcile.button("⟳ Reconcile")
auto_refresh = toggle_col.toggle("Auto-refresh (30 s)")

if do_fetch:
    backend, err = _try_backend()
    if err:
        st.error(f"Backend not configured: {err}")
    else:
        args = types.SimpleNamespace(job_id=None)
        buf = io.StringIO()
        conn = _open_db()
        try:
            with st.spinner("Fetching…"):
                run_fetch(args, settings, backend, conn, stdout=buf)
        except Exception as exc:
            st.error(str(exc))
        else:
            st.toast(buf.getvalue().strip() or "fetch done")
        finally:
            conn.close()
        st.rerun()

if do_reconcile:
    backend, err = _try_backend()
    if err:
        st.error(f"Backend not configured: {err}")
    else:
        args = types.SimpleNamespace()
        buf = io.StringIO()
        conn = _open_db()
        try:
            with st.spinner("Reconciling…"):
                run_reconcile(args, settings, backend, conn, stdout=buf)
        except Exception as exc:
            st.error(str(exc))
        else:
            st.toast(buf.getvalue().strip() or "reconcile done")
        finally:
            conn.close()
        st.rerun()

# ── metrics ───────────────────────────────────────────────────────────────────
conn = _open_db()
try:

    def _count(status: str) -> int:
        return conn.execute("SELECT count(*) FROM jobs WHERE status=?", (status,)).fetchone()[0]

    n_pending = _count("pending")
    n_running = _count("running")
    n_complete = _count("complete")
    n_error = _count("error")

    job_rows = conn.execute(
        "SELECT id, smiles, method, basis, job_kind, status, "
        "submitted_at, started_at, completed_at, ssh_jobid, error_msg "
        "FROM jobs ORDER BY id DESC LIMIT 200"
    ).fetchall()

    wf_rows = conn.execute(
        "SELECT id, kind, smiles, status, created_at, summary_json "
        "FROM workflows ORDER BY id DESC LIMIT 50"
    ).fetchall()
finally:
    conn.close()

m1, m2, m3, m4 = st.columns(4)
m1.metric("Pending", n_pending)
m2.metric("Running", n_running, delta=None)
m3.metric("Complete", n_complete)
m4.metric("Error", n_error)

st.divider()

# ── job table ─────────────────────────────────────────────────────────────────
_STATUS_CSS = {
    "complete": "color: #28a745; font-weight: bold",
    "running": "color: #fd7e14; font-weight: bold",
    "pending": "color: #007bff; font-weight: bold",
    "error": "color: #dc3545; font-weight: bold",
}


def _status_style(val: str) -> str:
    return _STATUS_CSS.get(val, "")


if not job_rows:
    st.info("No jobs in the database.")
else:
    df = pd.DataFrame([dict(r) for r in job_rows])
    # Truncate long SMILES so the table stays readable.
    df["smiles"] = df["smiles"].apply(lambda s: s if len(s) <= 30 else s[:27] + "…")
    st.dataframe(
        df.style.map(_status_style, subset=["status"]),
        use_container_width=True,
        hide_index=True,
    )

# ── workflows ─────────────────────────────────────────────────────────────────
if wf_rows:
    with st.expander(f"Workflows ({len(wf_rows)})"):
        wf_df = pd.DataFrame([dict(r) for r in wf_rows])
        st.dataframe(
            wf_df.style.map(_status_style, subset=["status"]),
            use_container_width=True,
            hide_index=True,
        )

# ── auto-refresh countdown ───────────────────────────────────────────────────
if auto_refresh:
    placeholder = st.empty()
    for remaining in range(30, 0, -1):
        placeholder.caption(f"Next refresh in {remaining} s…  (toggle off to stop)")
        time.sleep(1)
    placeholder.empty()
    st.rerun()
