# oscpipe

[![DOI](https://zenodo.org/badge/1246342088.svg)](https://zenodo.org/badge/latestdoi/1246342088) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

An end-to-end DFT pipeline for organic-semiconductor (OSC) molecules: from a SMILES
string to an optimized geometry and electronic properties, with the heavy Gaussian 16
jobs dispatched to a remote HPC cluster over SSH. Ships a CLI and a Streamlit dashboard,
backed by a SQLite job/result store.

## Features

- **SMILES → properties** in one command: RDKit 3D embed → Gaussian 16 input → optimize → HOMO / LUMO / gap.
- **Workflows:** ground-state properties; reorganization energy λ_h (4-point) + Marcus hopping rate; TD-DFT / UV–Vis excited states; oligomer-length sweep with 1/n extrapolation to the polymer limit; batch screening.
- **Remote dispatch:** key-auth SSH to a UGE/`qsub` cluster (Gaussian runs there); a local backend for development.
- **Resumable:** every job and result is recorded in SQLite; `fetch` / `reconcile` recover asynchronous runs.
- **Dashboard:** a Streamlit app with a properties view and a live workstation status board.

## Requirements

- Python ≥ 3.11
- **Gaussian 16 (`g16`)** on a remote host reachable via SSH key auth — Gaussian is commercial and is **not** bundled. The local dev backend runs without it (no real DFT).
- Host and scratch directory are set in `config.local.toml` (gitignored; see `docs/RUNBOOK.md`).

## Quick start

```bash
python3 -m venv .venv
.venv/bin/pip install -e ".[dev]"
oscpipe preflight                 # checks g16, qstat, scratch dir, DB
oscpipe submit "c1ccsc1"          # thiophene → Properties workflow
oscpipe fetch                     # pull completed logs back + parse
oscpipe status                    # print the job/result table
```

Example `status` output (values illustrative):

```
  id  status      method/basis           chg      HOMO     gap  smiles
   1  complete    b3lyp/6-31g*             0    -6.210   5.740  c1ccsc1
```

## Commands

| Command | Purpose |
|---|---|
| `submit <smiles>` | One molecule → ground-state properties (HOMO / LUMO / gap) |
| `batch <file>` | Screen a list (CSV with a `smiles` column, or `-` for stdin) |
| `lambda <smiles>` | Reorganization energy λ_h (4-point) + Marcus rate |
| `uvvis <smiles>` | TD-DFT excited states (UV–Vis); `--from-log` reuses an optimized geometry |
| `oligomer <repeat-unit>` | n-mer sweep + 1/n extrapolation to the polymer limit |
| `status` / `fetch` / `reconcile` | Job table / pull results / sync with the queue |
| `preflight` | Remote-host health check (g16, qstat, scratch, DB) |

`submit`, `batch`, `uvvis`, `oligomer` accept `--method` / `--basis` (defaults `b3lyp` / `6-31g*`).

## Documentation

- `docs/ARCHITECTURE.md` — layered design (chem · dft · dispatch · store · workflows · cli/app)
- `docs/RUNBOOK.md` — day-to-day usage and host configuration
- `docs/adr/` — architecture decision records

## Tests

```bash
.venv/bin/python -m pytest        # 133 passed, 2 skipped (real-workstation markers)
```

## License

MIT — see [`LICENSE`](LICENSE).
