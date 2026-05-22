# oscpipe

OSC DFT pipeline for the Smart Ring Prototype. Replaces `OSC-pipeline/` once the
cutover checklist passes (see `docs/CUTOVER_CHECKLIST.md`).

## Quick start

```bash
python3 -m venv .venv
.venv/bin/pip install -e ".[dev]"
oscpipe preflight
oscpipe submit "c1ccsc1"
```

See `docs/ARCHITECTURE.md` for the layered design and `docs/RUNBOOK.md` for
day-to-day usage.

## Status

All five Smart Ring workflows (WF1–WF5) and the reconciliation gate (R1) pass
end-to-end against the workstation `user@203.0.113.10`. See
`docs/CUTOVER_CHECKLIST.md` for verification details. Tests: 106 pass, 2 skipped
(real-workstation markers).
