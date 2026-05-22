# oscpipe — Runbook

## First-time setup (laptop)

1. Create venv and install:
   ```bash
   python3 -m venv .venv
   .venv/bin/pip install -e ".[dev]"
   ```
2. Pin the workstation host key (one-time):
   ```bash
   ssh user@203.0.113.10     # accept fingerprint, then close
   cp ~/.ssh/known_hosts ./known_hosts.local
   ```
3. Install your SSH public key on the workstation:
   ```bash
   ssh-copy-id user@203.0.113.10
   ```
4. Create `config.local.toml` next to where you run the CLI:
   ```toml
   remote_host       = "203.0.113.10"
   remote_user       = "user"
   remote_key_file   = "~/.ssh/id_ed25519"
   known_hosts_path  = "./known_hosts.local"
   remote_work_dir   = "/home/user/oscpipe_work"
   ```
5. Verify:
   ```bash
   oscpipe preflight
   ```

## Daily use

- Submit one: `oscpipe submit "c1ccsc1" --method b3lyp --basis 6-31g*`
- Batch: `oscpipe batch smiles.csv`
- λ_reorg + Marcus: `oscpipe lambda "c1ccsc1"`
- UV-Vis: `oscpipe uvvis "c1ccsc1" --nstates 10`
- Streamlit GUI: `streamlit run app/Home.py`

## When something breaks

- Job stuck `pending` → `oscpipe reconcile` syncs with qstat
- Workstation dropped during a batch → reconcile recovers orphans; re-submit
  the still-missing SMILES
- `.log` parses with NaN values → check `is_log_complete()`; rerun with the
  same submit (cache lookup is by signature, not log file)
