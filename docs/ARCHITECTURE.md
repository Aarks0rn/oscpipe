# oscpipe — Architecture

Single-purpose pipeline for the Smart Ring Prototype. Replaces OSC-pipeline/.

## Layers

```
chem      → SMILES sanitisation + 3D embedding (rdkit primary, ob fallback)
dft       → Gaussian .com writer + .log parser (pure functions, no IO orchestration)
dispatch  → Backend protocol; SSH-only for first cut (paramiko, key-only, strict known_hosts)
store     → SQLite, two tables (jobs + results) + workflows + view, explicit schema.sql
analysis  → Marcus (λ_reorg + rate), Indo (transfer integral), UV-Vis broadening
cli       → `oscpipe submit | batch | lambda | uvvis | status | fetch | reconcile | preflight`
app/      → Streamlit multipage: 1_Properties.py + 2_Workstation.py
```

No imports across layers go upward. `cli` and `app` are the only entry points.

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
