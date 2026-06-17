# ADR 0006 — Colocated/local execution reuses `SshBackend` over loopback, not a native async local backend

**Status:** Accepted
**Date:** 2026-06-17
**Authors:** Aarks0rn (improve-codebase-architecture session)

## Context

The laptop runs the `oscpipe` driver and polls `g16` on the workstation
(`MiniLab`) over a long-lived paramiko connection. A multi-hour sweep dropped it:
a diF oligomer sweep crashed mid-poll with `SSHException("SSH session not active")`
(job 178), leaving its n=2 opt as an orphan `g16`. Two problems surfaced:

1. **The poll connection is a single point of failure.** `_run`'s reconnect only
   caught `(OSError, EOFError)`, not paramiko's `SSHException` — a transport that
   has gone inactive crashed the whole sweep (fixed separately: `SSHException`
   added to the retry tuple in `dispatch/ssh.py`).
2. **The driver is tied to the laptop.** Even with the reconnect fix, a laptop
   sleep / network move kills the driver. A multi-day DFT campaign should not
   depend on the laptop staying awake and connected.

The fix for (2) is to **colocate** the driver on the workstation — run `oscpipe`
on `MiniLab` so the poll loop is local. That raised the question of *which
dispatch backend* a colocated driver uses.

The dispatch seam (`_make_backend`, selected by `settings.backend`) has two
adapters whose **interfaces diverge** (not just method names — the async-ness,
detach, concurrency, and reattach invariants the `run` layer depends on):

| contract | `SshBackend` | `LocalBackend` |
|---|---|---|
| `submit` | async — `nohup` detach, returns a **live PID** | synchronous — blocks until `g16` exits, returns the label |
| survives driver death | yes (`nohup`) | no (inline `subprocess.run`) |
| `poll` | liveness (`kill -0`) + log | log only (the job already finished) |
| concurrency (`resolve_concurrent`, `max_lanes`) | yes | no (submit blocks) |
| reattach (ADR-0005) | yes (PID + signature) | n/a |
| `cancel` | `kill` | no-op |

`run`-layer machinery — `resolve_concurrent` lane parallelism and
`_reattach_or_submit` (ADR-0005) — assumes the **async** contract. `LocalBackend`
is a synchronous dev/test double (its own docstring: "For tests + dev only") and
cannot honour it. So a colocated driver needs a *local* path with `SshBackend`'s
async/detach/PID/reattach contract.

## Decision

**Colocated/local execution uses `SshBackend` over loopback** —
`backend = "ssh"`, `remote_host = "localhost"` in the workstation-side config.
A loopback SSH connection reuses the full, proven async contract
(detach / live-PID poll / `resolve_concurrent` lanes / `_reattach_or_submit`)
for **zero new code**, and a loopback socket does not drop the way a laptop→host
network link does (the `SSHException` reconnect from the fix above covers any
residual edge).

**Do not build a native async `LocalDirectBackend`** (subprocess `nohup` +
local-PID poll + `_jobs` map + reattach hooks). It would re-implement, line for
line, the contract `SshBackend` already ships.

`LocalBackend` **stays a synchronous dev/test double** — it is not promoted to a
production/colocation backend.

## Consequences

- **Cutover config (workstation):** `backend = "ssh"`, `remote_host = "localhost"`;
  `MiniLab` SSHes to itself (gen a key, add its public half to `MiniLab`'s
  `~/.ssh/authorized_keys`). All write paths — `db_path` / `work_dir` /
  `log_dir` / `scratch_dir` (`GAUSS_SCRDIR`) — and the repo + venv live under
  `/home/MiniLab/Lab-Data/SA26`, so the colocated driver never writes outside
  that scope.
- The network-SSH failure class (a laptop-side drop crashing the driver) is
  removed at the root, not just recovered from. The driver no longer depends on
  the laptop being awake/connected.
- No new dispatch adapter, no new test surface. The seam keeps two adapters, but
  the production/colocation path is `SshBackend` (one async contract over two
  transports: network and loopback); `LocalBackend`'s divergent synchronous
  semantics stay confined to dev/tests.
- **Known gap (accepted):** the testability friction that motivated the candidate
  remains — tests on the synchronous `LocalBackend` still do not exercise the
  async poll/reattach/lane contract prod uses. We accept it rather than pay for a
  third adapter; the async contract is covered by the `SshBackend` unit tests
  (incl. the `SSHException` reconnect test) and by live runs.

## Alternatives rejected

- **Native async `LocalDirectBackend`** (subprocess instead of paramiko).
  ~50–60 lines plus a test surface that duplicates `SshBackend`'s detach / PID
  poll / reattach. Rejected — loopback reuses that contract for free; the only
  thing it buys is "no SSH dependency at all," which is not worth the duplication.
- **Deepen `LocalBackend` in place to the async contract.** Breaks its sole role
  as the synchronous test double the suite relies on. Rejected.
- **Keep the laptop driver over network SSH and rely only on the `SSHException`
  reconnect.** Hardens the drop but does not decouple the campaign from laptop
  connectivity (sleep / network move still kills the driver). Colocation removes
  that dependency; the reconnect fix stays as defence-in-depth for the loopback
  path.

## References

- `src/oscpipe/dispatch/ssh.py` — `SshBackend`; `_run` reconnect (now catches
  `SSHException`)
- `src/oscpipe/dispatch/local.py` — `LocalBackend` (synchronous dev/test double)
- `src/oscpipe/cli/main.py` — `_make_backend` (selects on `settings.backend`)
- `src/oscpipe/run.py` — `resolve_concurrent`, `_reattach_or_submit`
- ADR 0005 — orphan-prevention / reattach (the async contract this relies on)
- ADR 0003 — the `JobRunner.resolve` seam
- `settings.remote_host` — set to `localhost` for the colocated config
