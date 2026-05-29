# ADR 0003 — One job-resolution seam (`JobRunner.resolve`); persist-by-kind is total; λ_h reuses the cache

**Status:** Accepted
**Date:** 2026-05-29
**Authors:** Aarks0rn (improve-codebase-architecture session)

## Context

After ADR 0001 split orchestration into `run` + `workflows`, two parallel
implementations of "run one DFT job end-to-end" remained:

1. **CLI single-shot** (`cli.run_submit`, `cli.run_uvvis`): canonicalise → signature
   → cache-check → write `.com` → `JobRunner.submit_poll_once` → (async? return) →
   `finish_completed`. The cache-check + submit + error/finish tail was duplicated
   almost verbatim across the two commands.
2. **The λ_h workflow** (`workflows.lambda_h`): no cache-check →
   `JobRunner.submit_and_wait` → parse energies **inline** with
   `gaussian.parse_properties`.

The two paths diverged on three axes — caching, sync model (poll-once vs block),
and persistence — and that divergence caused real friction:

- **No single home for "job_kind → parse + persist."** `finish_completed`
  dispatched only `tddft` vs *everything-else-as-properties* (2 of the schema's 6
  kinds); the other four were parsed inline in `lambda_h`. A λ_h sub-job that was
  detached (Ctrl-C) and later resumed through `fetch`/`reconcile` was therefore
  mis-persisted — a `sp_dimer` log was parsed as monomer *properties*, writing a
  bogus HOMO/LUMO row, and the workflow never produced J. The detach message even
  told the user to run `fetch`, which could not actually finish the workflow.
- **λ_h never consulted the signature cache.** It computed each sub-job's signature
  (and stored it on the row) but never called `find_complete_by_signature`, so
  re-running `oscpipe lambda` after any failure re-submitted all five multi-hour
  DFT jobs even when identical completed jobs existed.

## Decision

Introduce **`JobRunner.resolve(job, label, build_com, *, wait, need_log=False)`** —
the single operation "resolve one job to its result":

```
signature cache hit  → return stored result (build_com is never called, so the
                       expensive embed/.com write is skipped)
else                 → submit → poll once (wait=False) | block (wait=True)
                       → fetch → parse → persist by job_kind → return result
```

- **persist-by-kind (`_persist_by_kind`) is the one home** for the job_kind → parse
  + persist mapping, and is **total** over the kinds the schema defines:
  `properties` / `sp_neutral` / `sp_cation` → a results row; `tddft` → a spectra
  row; `sp_dimer` → status + log_path only (its J_hole is a workflow-level quantity,
  derived by `analysis.indo` into `workflows.summary_json`); an unhandled kind
  (e.g. `freq`, never produced yet) raises rather than persisting silently.
- `cli.run_submit` / `cli.run_uvvis` become thin `resolve` adapters sharing one
  exit-code tail (`_submit_rc`).
- `workflows.lambda_h` calls `resolve` for all five jobs — gaining the cache — and
  reads energies from the returned result instead of re-parsing the logs.
- **`build_com` is a thunk** so the embed/write is skipped on a cache hit.
- **`need_log=True`** invalidates a cache hit whose local `.log` is gone (the caller
  needs the geometry or the dimer eigenvalues), forcing a recompute.

**Resume scope.** Re-running `oscpipe lambda` is the resume path: cached sub-jobs are
skipped and only the missing ones run. `fetch`/`reconcile` persist a resumed
sub-job's *raw* result correctly (via the same `_persist_by_kind`) but do **not**
recompute workflow-level quantities — that would need a workflow state machine /
checkpointing, a feature beyond this deepening.

## Consequences

- The single-shot submit arc is no longer duplicated; the job_kind → persist mapping
  has one home. A resumed `sp_dimer` no longer writes a bogus properties row
  (`tests/test_run.py::test_fetch_resumes_dimer_without_bogus_results_row`).
- **DB content change (intentional, approved):** λ_h sub-jobs now get results rows —
  `sp_neutral` / `sp_cation` carry `energy_ev`; `sp_dimer` still has none. The two
  opt sub-jobs already carried `job_kind='properties'`, so populating their results
  made them surface in the Properties dashboard (cation_opt would show its +1 alpha
  orbitals under HOMO/LUMO with no charge cue). The Properties browser query now
  excludes `lambda_h` workflow sub-jobs (`app/pages/1_Properties.py`); they remain
  on the Workstation page and in the workflow summary. The SP/dimer jobs were
  already filtered out by `job_kind`. `oscpipe status` keeps every job (it is the
  raw monitor) but gained a `chg` column so a charged sub-job's HOMO/gap reads as a
  cation result, not a neutral one. **Rule:** the Properties *results browser* hides
  λ_h sub-jobs; `status` and the Workstation page show every job, with charge /
  job_kind for context.
- λ_h re-run reuses the cache: 5 jobs on the first run, 0 new on re-run, and the
  re-run reproduces the same λ_h (`test_lambda_rerun_reuses_cache_no_new_jobs`).
- A cached re-run inserts a fresh `lambda_h` workflow row with a correct
  `summary_json`, but its five sub-jobs stay linked to the *first* run's
  `workflow_id` (a cache hit neither inserts nor relinks a job row). The second
  workflow thus has a correct summary but no member jobs — accepted as the cost of
  "re-run = cheap reproduce", not a state machine.
- `submit_poll_once`, `submit_and_wait`, and `RunResult` are removed — subsumed by
  `resolve` / `Resolved`.
- 112 prior tests stay green, 5 added (117 total / 2 skipped). λ_h is bit-identical
  on the known case — the 2.72114 eV regression test is untouched.
- **Minor, deliberate:** a parse-failure-after-complete now exits 1 (was an
  undocumented, untested 2); the cache-hit log lines were simplified (only the
  `"cache hit"` substring was ever contracted by a test).

## Alternatives rejected

- **Keep the cache-check in each caller; `resolve` does only submit + persist.**
  Re-duplicates the cache-check into `lambda_h` — exactly what this removes. Rejected.
- **Pass `com_text` eagerly instead of a `build_com` thunk.** Regresses the
  embed-skip-on-cache-hit optimisation, which matters for batch / screening re-runs.
  Rejected.
- **Make `fetch`/`reconcile` finish a λ_h workflow** (recompute λ_h/J when the last
  sub-job lands). Needs a workflow state machine + checkpointing — a feature, not a
  deepening; re-running `lambda` over the cache achieves resume cheaply. Rejected (scope).
- **Give `sp_dimer` a results row too** (full totality). The schema has no
  J / HOMO-1 column and J is a workflow-level quantity; a per-job dimer row would be
  misleading. Rejected.

## References

- `docs/ARCHITECTURE.md` — `run` layer (`JobRunner.resolve`)
- `src/oscpipe/run.py` — `resolve`, `_persist_by_kind`, `_load_cached`, `Resolved`
- `src/oscpipe/cli/main.py` — `run_submit` / `run_uvvis` adapters + `_submit_rc`
- `src/oscpipe/workflows/lambda_h.py` — five `resolve` calls; energies from results
- `app/pages/1_Properties.py` — Properties-browser filter (excludes λ_h sub-jobs)
- `tests/test_run.py` — cache hit, `need_log` recompute, total persist, λ_h re-run,
  Properties-browser filter
- ADR 0001 — the `run` / `workflows` split this builds on
