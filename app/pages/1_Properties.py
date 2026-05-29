"""Properties page — single SMILES + batch submit, results browser + UV-Vis."""

from __future__ import annotations

import io
import json
import types

import numpy as np
import pandas as pd
import streamlit as st

from oscpipe.cli.main import _make_backend, run_batch, run_submit
from oscpipe.settings import load
from oscpipe.store import db

st.title("Properties")
st.caption("HOMO / LUMO / gap / dipole")

settings = load()


def _open_db():
    return db.open(settings.db_path)


def _try_backend():
    try:
        return _make_backend(settings), None
    except ValueError as e:
        return None, str(e)


# ── single SMILES submit ──────────────────────────────────────────────────────
with st.expander("Submit single SMILES", expanded=True):
    with st.form("single"):
        smiles_in = st.text_input("SMILES", placeholder="c1ccsc1")
        c1, c2, c3, c4 = st.columns(4)
        method = c1.text_input("Method", value="b3lyp")
        basis = c2.text_input("Basis", value="6-31g*")
        charge = c3.number_input("Charge", value=0, step=1)
        mult = c4.number_input("Mult", value=1, min_value=1, step=1)
        go = st.form_submit_button("Submit")

    if go and smiles_in.strip():
        backend, err = _try_backend()
        if err:
            st.error(f"Backend not configured: {err}")
        else:
            args = types.SimpleNamespace(
                smiles=smiles_in.strip(),
                method=method,
                basis=basis,
                charge=int(charge),
                mult=int(mult),
            )
            buf = io.StringIO()
            conn = _open_db()
            try:
                with st.spinner("Submitting…"):
                    rc = run_submit(args, settings, backend, conn, stdout=buf)
            except Exception as exc:
                st.error(str(exc))
            else:
                msg = buf.getvalue() or ("submitted" if rc == 0 else "error")
                (st.success if rc == 0 else st.error)(msg)
            finally:
                conn.close()

# ── batch submit ─────────────────────────────────────────────────────────────
with st.expander("Batch submit (one SMILES per line)"):
    batch_raw = st.text_area("SMILES list", height=110, placeholder="c1ccsc1\nc1ccccc1")
    bc1, bc2 = st.columns(2)
    b_method = bc1.text_input("Method##batch", value="b3lyp")
    b_basis = bc2.text_input("Basis##batch", value="6-31g*")

    if st.button("Submit batch"):
        lines = [ln.strip() for ln in batch_raw.splitlines() if ln.strip()]
        if not lines:
            st.warning("Enter at least one SMILES.")
        else:
            backend, err = _try_backend()
            if err:
                st.error(f"Backend not configured: {err}")
            else:
                args = types.SimpleNamespace(
                    file="-", method=b_method, basis=b_basis, charge=0, mult=1
                )
                buf = io.StringIO()
                conn = _open_db()
                try:
                    with st.spinner(f"Submitting {len(lines)} molecules…"):
                        rc = run_batch(
                            args,
                            settings,
                            backend,
                            conn,
                            stdin=io.StringIO("\n".join(lines)),
                            stdout=buf,
                        )
                except Exception as exc:
                    st.error(str(exc))
                else:
                    msg = buf.getvalue()
                    (st.success if rc == 0 else st.warning)(msg)
                finally:
                    conn.close()

# ── results table ─────────────────────────────────────────────────────────────
st.divider()

hdr_col, btn_col = st.columns([6, 1])
hdr_col.subheader("Results")
if btn_col.button("↺ Refresh"):
    st.rerun()

conn = _open_db()
try:
    # Exclude λ_h workflow sub-jobs: neutral_opt / cation_opt carry
    # job_kind='properties' but are workflow internals (cation_opt would show its
    # +1 alpha orbitals under HOMO/LUMO with no charge cue). Batch jobs also carry
    # a workflow_id but are genuine standalone requests, so filter on kind only.
    prop_rows = conn.execute(
        "SELECT id, smiles, method, basis, status, submitted_at, "
        "homo_ev, lumo_ev, gap_ev, dipole_debye "
        "FROM v_jobs_with_results WHERE job_kind = 'properties' "
        "AND (workflow_id IS NULL OR workflow_id NOT IN "
        "(SELECT id FROM workflows WHERE kind = 'lambda_h')) "
        "ORDER BY id DESC"
    ).fetchall()
    tddft_rows = conn.execute(
        "SELECT j.id, j.smiles, j.status, r.spectra_json "
        "FROM jobs j LEFT JOIN results r ON r.job_id = j.id "
        "WHERE j.job_kind = 'tddft' AND j.status = 'complete' ORDER BY j.id DESC"
    ).fetchall()
finally:
    conn.close()

_STATUS_CSS = {
    "complete": "color: #28a745; font-weight: bold",
    "running": "color: #fd7e14; font-weight: bold",
    "pending": "color: #007bff; font-weight: bold",
    "error": "color: #dc3545; font-weight: bold",
}


def _status_style(val: str) -> str:
    return _STATUS_CSS.get(val, "")


if not prop_rows:
    st.info("No property jobs yet. Submit a SMILES above.")
else:
    df = pd.DataFrame([dict(r) for r in prop_rows])
    for col in ("homo_ev", "lumo_ev", "gap_ev", "dipole_debye"):
        df[col] = df[col].apply(lambda v: f"{v:.3f}" if v is not None else "—")
    st.dataframe(
        df.style.map(_status_style, subset=["status"]),
        use_container_width=True,
        hide_index=True,
    )

# ── UV-Vis spectra ────────────────────────────────────────────────────────────
if tddft_rows:
    st.divider()
    st.subheader("UV-Vis (TDDFT)")

    options = {f"job {r['id']} — {r['smiles']}": r for r in tddft_rows}
    sel_key = st.selectbox("Select job", list(options.keys()))
    sel = options[sel_key]

    if sel["spectra_json"]:
        states = json.loads(sel["spectra_json"])

        # Gaussian broadening of the stick spectrum.
        wl = np.linspace(200, 800, 601)
        sigma = 20.0  # nm
        spectrum = np.zeros_like(wl)
        for s in states:
            spectrum += s["f"] * np.exp(-0.5 * ((wl - s["wavelength_nm"]) / sigma) ** 2)

        chart_df = pd.DataFrame({"Wavelength (nm)": wl, "Intensity": spectrum})
        st.line_chart(chart_df.set_index("Wavelength (nm)"), height=260)

        sticks = pd.DataFrame(states).rename(
            columns={
                "n": "state",
                "wavelength_nm": "λ (nm)",
                "energy_ev": "E (eV)",
                "f": "osc. strength",
            }
        )
        st.dataframe(
            sticks[["state", "λ (nm)", "E (eV)", "osc. strength"]],
            use_container_width=True,
            hide_index=True,
        )
