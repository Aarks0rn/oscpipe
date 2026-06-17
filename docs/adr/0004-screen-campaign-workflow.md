# ADR 0004 — `oscpipe screen`: campaign as sequential composition over the signature cache; store read API + `oscpipe export`

**Status:** Accepted
**Date:** 2026-06-11
**Authors:** Aarks0rn (improve-codebase-architecture session)

## Context

The roadmap's M3 screen needs 5–10 candidates run through the full pipeline
(properties → TD-DFT → λ_h) over ~2–3 months of workstation time. Today that
is manual: per candidate, the operator runs `submit`, `uvvis`, `lambda` by
hand and tracks progress in their head or in shell scripts
(`polymer-donor-db/scripts/submit-*.sh`). Resume after an interruption
(workstation reset, SSH death, Ctrl-C) depends on the operator knowing what
already ran. ADR 0003 made every sub-job cache-resumable, but nothing
composes the three stages per candidate.

On the read side, results leave the DB via ad-hoc raw SQL in project scripts
— the schema is an implicit interface nobody owns — and the M1 dataset
release needs a flat results dump that does not exist.

The ARCHITECTURE.md scope gate requires an ADR for features beyond the four
core workflows. Screen and export compose/reshape those workflows' outputs;
neither adds physics.

## Decision

- **`workflows/screen.py::run_screen`** — per candidate, sequential and
  blocking: properties (`resolve wait=True`) → TD-DFT (`resolve wait=True`) →
  `run_lambda_h` (unmodified). Any failure records the candidate and
  continues with the next one. One `workflows` row `kind='screen'` whose
  `summary_json` carries per-candidate progress
  (`properties`/`tddft`/`lambda_h` per entry); child λ_h workflows keep their
  own rows and jobs stay linked to them.
- **Resume model: re-run the same CSV.** Completed sub-jobs are signature
  cache hits (`"cached"` in the new summary); only missing work submits.
  No state machine, no checkpoint table. A re-run inserts a fresh screen row
  — same accepted cost as ADR 0003's cached re-run.
- **Signature alignment is deliberate.** Screen's properties step defaults to
  b3lyp/6-31g\*\* so λ_h's `neutral_opt` signature equals it — λ_h cache-hits
  screen's own opt (one multi-hour job saved per candidate). The TD-DFT step
  byte-matches `oscpipe uvvis` signatures, so prior uvvis runs cross-hit.
  Consequence: screen's basis default (6-31g\*\*) differs from
  `submit`/`batch`'s 6-31g\* — a prior batch at defaults will NOT cross-hit.
- **`run.make_job`** is the single home of the canonical-SMILES → signature →
  label → pending-Job convention; `run_submit`, `run_uvvis` and screen use it.
  `lambda_h`/`oligomer_sweep` deliberately do not: their signature job_kind
  differs from the row job_kind (ADR 0003) and routing them through the
  helper would silently change signatures and orphan the live cache
  (guard test: `test_make_job_matches_existing_signature_convention`).
- **Read side:** `store.db` gains `list_complete_oligomer_sweeps`,
  `candidate_summary`, `list_candidate_smiles`; `oscpipe export --csv` is a
  thin adapter over them (one flat CSV, no other formats). External scripts
  consume the helpers instead of raw SQL.

## Consequences

- One command (`oscpipe screen candidates.csv`) covers the M3 campaign; after
  any interruption the same command resumes it. Progress is queryable from
  the screen row's summary instead of reconstructed from job rows.
- A detached candidate (Ctrl-C mid-wait) is recorded with its raw status and
  counts as failed for that run; the re-run picks it up via the cache.
- Screen's properties sub-jobs appear in the Properties browser (they are
  genuine neutral results); λ sub-jobs stay hidden because they link to the
  child lambda_h workflow ids. No app changes needed.
- `oscpipe export` gives M1 a flat one-row-per-candidate CSV;
  `protocol_offset.py` now reads sweeps through the store API, so the
  workflows-table shape is encoded in one place.

## Alternatives rejected

- **Workflow state machine / checkpointing.** Re-rejected (ADR 0003);
  cache-resume is sufficient and free.
- **Parallel fire-and-forget submission per candidate.** λ_h is internally
  sequential (each SP needs the opposite opt's geometry); blocking keeps
  failure handling and queue load simple; serial restart is cheap via cache.
- **TD-DFT at the optimised geometry from the properties step.** Would
  decouple screen's TD-DFT from `oscpipe uvvis` cache entries and require the
  geometry discriminator on every screen run. Planar-embed TD-DFT matches
  current uvvis semantics. Revisit if accuracy demands.
- **Screen inside `cli` like `batch`.** It is a multi-job scientific
  sequence; ADR 0001 puts those in `workflows`.
- **New campaign/candidate tables; JSON/report export formats.** YAGNI —
  `workflows.summary_json` and one CSV suffice.

## References

- `src/oscpipe/workflows/screen.py` — the workflow
- `src/oscpipe/run.py` — `make_job`
- `src/oscpipe/store/db.py` — read helpers; `src/oscpipe/cli/main.py` — `run_screen`, `run_export`
- `tests/test_cli_screen.py`, `tests/test_run.py::test_make_job_matches_existing_signature_convention`
- ADR 0001 (layers), ADR 0003 (resolve seam, cache-resume)
- Workspace `ROADMAP.md` — M1 dataset, M3 screen
