# ADR 0005 — Orphan UGE jobs: prevent by reattach, reclaim by a read-only lister (never auto-qdel)

**Status:** Accepted
**Date:** 2026-06-11
**Authors:** Aarks0rn (scrutinize session)

## Context

Re-running a sweep / λ_h / screen kept leaving orphan UGE jobs — remote jobs
still burning the workstation's slots that no local process tracks. Two
mechanisms, found by tracing `run.resolve` → `dispatch.ssh`:

1. **Detach → resubmit (dominant).** Ctrl-C / a laptop disconnect during
   `JobRunner._wait` returns `"detached"` (run.py — intentional: a multi-hour
   job should survive a disconnect). The job row stays `running` with its
   `ssh_jobid`, and the UGE job keeps running. On re-run, `resolve`'s cache
   check (`find_complete_by_signature`) matches only `status='complete'`, so the
   still-`running` row is *not* a hit → a **duplicate UGE job is submitted**,
   orphaning the first.
2. **Stale-log false-complete (secondary).** `SshBackend._poll_log` trusts a
   `<label>.log` containing "Normal termination". A leftover log from a prior
   run makes a freshly-qsub'd job read as complete immediately, while the real
   UGE job runs on, orphaned. (Observed: DB job 20 "complete" 30 s after submit
   while UGE 71229 ran for real.)

Constraint (the workstation account `MyLab` is *shared* across lab members):
the tooling must never delete another person's job, and UGE **reuses job IDs**,
so an `ssh_jobid` alone is not a safe identity.

## Decision

**Prevent, don't clean up.** The fix is at submit time, deletes nothing:

- `run.JobRunner.resolve` gains a reattach step (`_reattach_or_submit`):
  before submitting, `store.db.find_inflight_by_signature` looks for a
  `pending`/`running` row with the same signature and a recorded `ssh_jobid`.
  If found, poll **that** job instead of qsub-ing a duplicate (and register its
  remote log path so a fresh process can poll/fetch it). `build_com` is not
  called on this path. This kills mechanism (1) at the root with zero deletion —
  it only ever touches jobs *we* recorded.

**Reclaim the rest with a read-only lister**, never an auto-`qdel`:

- `oscpipe orphans` lists redundant queued jobs and prints the `qdel` line for
  the user to run. A job is flagged only when **both** hold: its UGE name
  (from `qstat -u <user> -r`, untruncated) starts with one of our label stems
  `<slug>_<16hex signature>` (a lab-mate's hand-named job never matches), **and**
  that signature is already `complete` in our DB (so the running job is
  genuinely redundant). Name + DB-membership together are the ownership
  boundary; the raw reused ID is never trusted on its own.

This covers mechanism (2) and the give-up case (detached, never resumed)
without the tool ever deleting anything.

## Consequences

- Re-running after a detach reattaches (`reattach: job=… already queued — not
  resubmitting`) instead of duplicating. `resolve`'s contract widens: a
  same-signature in-flight job is now resumed, not re-submitted.
- `oscpipe orphans` is a diagnostic command alongside `status` / `reconcile` /
  `preflight` (read-only, no ADR-gated workflow scope). Verified against the
  live workstation: it surfaced orphan 71229 and correctly skipped the live
  n3_opt job (71232) whose signature was not yet complete.
- Not addressed (deliberate, scoped out): clearing stale remote logs *before*
  submit (the prevention for mechanism 2) — bigger, riskier change; the lister
  surfaces those orphans instead. The reset-killed-job case is unchanged
  (`reconcile` still marks those lost; they burn no slots).

## Alternatives rejected

- **Auto-`qdel` of orphans.** Rejected — the shared account + UGE ID reuse make
  any automated delete a risk to other members' jobs. The user reclaims slots
  by eye from the lister.
- **Trap Ctrl-C to `qdel` the in-flight job (turn detach into cancel).**
  Contradicts the intentional "detach survives disconnect" design and would
  throw away hours of compute on an accidental Ctrl-C. Reattach keeps the
  survival behavior and removes the duplicate instead.
- **Identify our jobs by `ssh_jobid` alone.** Unsafe — IDs are reused; matching
  must also check the name follows our label convention.

## References

- `src/oscpipe/run.py` — `_reattach_or_submit`; `store/db.py` —
  `find_inflight_by_signature`
- `src/oscpipe/dispatch/ssh.py` — `list_user_jobs`, `_parse_qstat_r`;
  `src/oscpipe/cli/main.py` — `_orphan_candidates`, `run_orphans`
- `tests/test_run.py` (reattach), `tests/test_cli_orphans.py` (lister + parser)
- ADR 0003 — the `resolve` seam this extends
