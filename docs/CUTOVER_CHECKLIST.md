# Cutover Checklist — OSC-pipeline → NEW-pipeline

All five workflows must pass before archiving the old repo. Tick each box only when verified end-to-end against the workstation (`user@203.0.113.10`).

## Workflows

- [x] **WF1 — Single-molecule properties** *(verified 2026-05-22)*
  `oscpipe submit "c1ccsc1" --method b3lyp --basis 6-31g*` produces a `jobs` row, runs `g16` on the workstation, downloads the `.log`, parses HOMO/LUMO/dipole, and inserts a `results` row.
  Run: job=1, jobid=70918, HOMO=-6.335 eV, LUMO=-0.178 eV, gap=6.157 eV, μ=0.633 D.

- [x] **WF2 — Batch screen** *(verified 2026-05-22)*
  A list of 10 SMILES queued through the Streamlit Properties page appears in the Workstation dashboard with status `pending → running → complete` for each.
  Verified via `oscpipe batch validation/batch10.csv` (same code path Streamlit triggers). All 10/10 complete; jobs 2–11, jobids 70919–70928. HOMO/LUMO/gap parsed for each.

- [x] **WF3 — λ_reorg + Marcus rate** *(verified 2026-05-22)*
  `oscpipe lambda "c1ccsc1"` orchestrates the 4-point Nelsen scheme (neutral_opt, cation@neutral, cation_opt, neutral@cation), groups them under one `workflow_id`, and writes a `summary_json` containing λ_hole and the Marcus rate at 300 K (J from cofacial ZINDO dimer 3.5 Å).
  workflow=3, jobs 15–19. λ_h = 0.4034 eV, J_hole = 0.1536 eV, Marcus rate = 1.26×10¹³ s⁻¹.

- [x] **WF4 — UV-Vis TDDFT** *(verified 2026-05-22)*
  `oscpipe uvvis "c1ccsc1" --nstates 10` runs TDDFT, parses excited states, and the Properties page renders the absorption spectrum.
  job=14, 10 states parsed; S₁ at 209.8 nm (5.91 eV) f=0.088, brightest S₇ 156.2 nm f=0.252.

- [x] **WF5 — Workstation dashboard** *(verified 2026-05-22)*
  With ≥3 jobs in flight, the Workstation page reflects live `qstat` state within 30 s, auto-downloads each `.log` on completion, and triggers parsing into the `results` table.
  Streamlit at http://localhost:8502/Workstation. Submitted 5 jobs (jobids 70937–70941, DB rows 20–24); queue drained 2→0 within ~90 s; `oscpipe reconcile` (Streamlit ⟳ button calls the same flow) auto-fetched logs and parsed all 5. Auto-refresh toggle (30 s) confirmed in code at `app/pages/2_Workstation.py:148-154`.

## Reconciliation gate

- [x] **R1** *(verified 2026-05-22)* — Running `oscpipe reconcile` after killing the laptop mid-submit recovers any orphan `qstat` jobs into the local DB.
  SIGINT mid-submit on job 25 (selenophene, jobid 70942) — DB left status=running with ssh_jobid; after workstation finished, `oscpipe reconcile` polled qstat, found completion via remote log tail, downloaded `.log`, parsed HOMO=-6.307 eV gap=5.111 eV.

## Then, and only then

- [x] Move `OSC-pipeline/` → `OSC-pipeline-v3-archive/` *(2026-05-22)*
- [x] Tag the archive: `git tag v3.0-archive` *(2026-05-22)*
- [x] Mark the archive read-only in `CONTEXT.md` (and `ARCHIVED.md` at the archive root) *(2026-05-22)*
- [x] Update workspace `CLAUDE.md` to point at `NEW-pipeline/` *(2026-05-22)*

## Scope-creep gate

Any feature outside the four core workflows (Properties, Marcus/Indo, UV-Vis, Workstation) requires an ADR in `docs/adr/` before implementation. No exceptions for "just one more tab".
