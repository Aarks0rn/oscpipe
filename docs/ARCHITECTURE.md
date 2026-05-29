# oscpipe — Architecture

Single-purpose pipeline for the Smart Ring Prototype. Replaces OSC-pipeline/.

## Layers

```
chem      → SMILES sanitisation + 3D embedding (rdkit primary, ob fallback);
            geometry I/O + π-stacked dimer construction (geometry.build_pi_stack_dimer)
dft       → Gaussian .com writer + .log parser (pure functions, no IO orchestration)
dispatch  → Backend protocol; SSH-only for first cut (paramiko, key-only, strict known_hosts)
store     → SQLite, two tables (jobs + results) + workflows + view, explicit schema.sql
analysis  → Marcus (λ_reorg + rate), Indo (transfer integral), UV-Vis broadening
run       → JobRunner.resolve: resolve one Gaussian job to its result — signature
            cache hit, else submit → poll/wait → fetch → parse → persist by
            job_kind. Single entry for the CLI single-shot commands and every
            workflow step; owns the remote-label + total persist-by-kind conventions.
workflows → multi-job sequences. `lambda_h.run_lambda_h` composes the 5 jobs
            (4-point Nelsen + π-stacked dimer SP) via JobRunner, then derives the
            physics through the `analysis` functions (no inline arithmetic).
cli       → `oscpipe submit | batch | lambda | uvvis | status | fetch | reconcile | preflight`
app/      → Streamlit multipage: 1_Properties.py + 2_Workstation.py
```

No imports across layers go upward. `cli` and `app` are the only entry points.
`run` composes `store` + `dft` + `dispatch` so callers stop threading
`(conn, settings, backend)` and re-deriving the label/job_kind conventions.
Both the single-shot commands and the workflow steps run through
`JobRunner.resolve`, so the signature cache and the persist-by-kind mapping are
defined once instead of per caller (see ADR 0003).
`workflows` sits on `run` + `analysis`: a workflow sequences jobs and owns no
physics of its own — the λ_h / transfer-integral / Marcus arithmetic lives in
`analysis` and is called, not re-implemented (the schema's `workflows` table
mirrors this layer).

## Scope

Smart Ring Prototype workflows only:

1. Properties (HOMO / LUMO / gap / dipole)
2. λ_reorg (4-point Nelsen) + Marcus rate
3. Transfer integral (ZINDO dimer)
4. UV-Vis TDDFT

Out of scope: herringbone, NICS, oligomer scan. Refer to the OSC-pipeline-v3-archive
if any of those become relevant for the Smart Ring path.

## Job model

Fire-and-forget. `oscpipe submit` writes a `jobs` row with status='pending' and
SSH-submits the .com to the workstation, then returns. Status moves
pending → running → complete | error as the Workstation dashboard polls qstat
and downloads completed logs.

## Security

SSH backend uses key auth only (`settings.remote_key_file`) and
`paramiko.RejectPolicy` against a pinned known_hosts file. Password mode is
disallowed at the type level (no `remote_password` field exists in `Settings`).

## Scope-creep gate

Any feature outside the four workflows above requires an ADR in `docs/adr/`
before implementation.
