# ADR 0001 — Orchestration lives in two layers (`run`, `workflows`), not in the CLI

**Status:** Accepted
**Date:** 2026-05-29
**Authors:** Aarks0rn (improve-codebase-architecture session)

## Context

`cli/main.py` had grown to ~1000 lines. Two distinct concerns were fused into it:

1. **Running one Gaussian job** — submit → poll/wait → fetch → parse → persist,
   plus the remote-label naming convention and a per-`job_kind` persist dispatch.
   This logic was forked (`_submit_job_and_poll_once` vs `_submit_and_wait`) and
   re-derived at every call site (`submit`, `uvvis`, `fetch`, `reconcile`, `lambda`).
2. **The λ_h workflow** — the 4-point Nelsen sequence + π-stacked dimer SP. Its
   physics was written inline in `run_lambda`, duplicating the tested
   `analysis.marcus.lambda_hole_from_4_points` and `analysis.indo.transfer_integral`
   (which had **zero production callers** — the shipping copy was the untested one).

The original layer list (`chem → dft → dispatch → store → analysis → cli → app`)
had no home for orchestration, so it accreted in `cli`.

## Decision

Introduce two layers below `cli`:

- **`run`** (`oscpipe/run.py`) — `JobRunner`: the lifecycle of **one** job. Owns the
  label convention and `job_kind` persist dispatch.
- **`workflows`** (`oscpipe/workflows/`) — multi-job sequences. `lambda_h.run_lambda_h`
  composes `JobRunner` calls and derives physics **through the `analysis` functions**,
  never inline. Mirrors the `workflows` table already in the schema.

New layer order: `chem → dft → dispatch → store → analysis → run → workflows → cli → app`.
No upward imports; `cli`/`app` remain the only entry points.

## Consequences

- The λ_h / transfer-integral arithmetic has a single home (`analysis`); the
  previously-dead tested functions now ship. λ_h on the known case is bit-identical
  (test asserts 2.72114 eV, still green).
- `cli.run_lambda` is a 3-line adapter; `submit`/`uvvis`/`fetch`/`reconcile` construct
  a `JobRunner` instead of threading `(conn, settings, backend)` and re-deriving labels.
- 112 tests green before and after (behaviour-preserving refactor). `app/` pages were
  verified to depend only on preserved public symbols (`_make_backend`, `run_*`).
- Dimer construction (`build_pi_stack_dimer`) moved to `chem.geometry`, concentrating
  π-stack geometry knowledge where the two open dimer bugs (clash → +40 eV, no-dispersion
  SP method) can be fixed locally rather than inside a 200-line CLI function.

## Alternatives rejected

- **Keep orchestration in `cli`, extract only a pure `lambda_h` composer.** Smaller, but
  leaves the job-runner fork and the `(conn, settings, backend)` threading in place.
  Rejected — chosen scope was the full layer split (see session).
- **A single `workflows/` layer holding both the runner and the sequences.** A single
  job is not a "workflow" (the schema distinguishes `jobs` from `workflows`); collapsing
  them would blur that. Rejected.

## References

- `docs/ARCHITECTURE.md` — layer list + the `run` / `workflows` descriptions
- Schema: `store/schema.sql` (`jobs`, `workflows` tables)
- Known dimer bugs: workspace `CLAUDE.md` → "Known issues / open items"
