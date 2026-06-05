"""Streamlit hub page. Routes to Properties and Workstation."""

import streamlit as st

st.set_page_config(page_title="oscpipe", layout="wide")
st.title("oscpipe — OSC DFT pipeline")
st.caption("SMILES → Gaussian 16 → properties")

st.markdown(
    """
    Pick a page from the sidebar:

    - **Properties** — submit single SMILES or a batch, view HOMO/LUMO/dipole results.
    - **Workstation** — live job dashboard (qstat, log download).

    Other analyses (λ_reorg / Marcus, UV-Vis, herringbone, NICS) are CLI / notebook
    only. See `docs/RUNBOOK.md`.
    """
)
