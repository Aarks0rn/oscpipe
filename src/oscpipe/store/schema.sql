-- oscpipe DB schema. Explicit, raw sqlite3.
-- Two tables: jobs (lifecycle) + results (scientific outputs).
-- A results row only exists for jobs.status='complete'.

PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS jobs (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    signature       TEXT    NOT NULL,        -- hash(canonical_smiles+method+basis+charge+mult)
    smiles          TEXT    NOT NULL,        -- canonicalised SMILES
    method          TEXT    NOT NULL,        -- e.g. b3lyp
    basis           TEXT    NOT NULL,        -- e.g. 6-31g*
    charge          INTEGER NOT NULL DEFAULT 0,
    mult            INTEGER NOT NULL DEFAULT 1,
    job_kind        TEXT    NOT NULL,        -- properties | tddft | freq | sp_cation | sp_neutral | sp_dimer
    status          TEXT    NOT NULL,        -- pending | running | complete | error | cancelled
    submitted_at    TEXT    NOT NULL,        -- ISO 8601 (local)
    started_at      TEXT,
    completed_at    TEXT,
    log_path        TEXT,                    -- local path to downloaded .log
    remote_log_path TEXT,                    -- workstation path
    ssh_jobid       TEXT,                    -- UGE qsub job id, if applicable
    error_msg       TEXT,
    workflow_id     INTEGER,                 -- groups multi-step workflows (λ_reorg = 4 jobs)
    notes           TEXT
);

CREATE INDEX IF NOT EXISTS idx_jobs_signature ON jobs(signature);
CREATE INDEX IF NOT EXISTS idx_jobs_status    ON jobs(status);
CREATE INDEX IF NOT EXISTS idx_jobs_workflow  ON jobs(workflow_id);

CREATE TABLE IF NOT EXISTS results (
    job_id        INTEGER PRIMARY KEY REFERENCES jobs(id) ON DELETE CASCADE,
    homo_ev       REAL,
    lumo_ev       REAL,
    gap_ev        REAL,
    dipole_debye  REAL,
    energy_ev     REAL,
    atoms_xyz     TEXT,                      -- XYZ-format string of optimised geometry
    spectra_json  TEXT                       -- TDDFT excited states / IR / etc, JSON blob
);

CREATE TABLE IF NOT EXISTS workflows (
    id            INTEGER PRIMARY KEY AUTOINCREMENT,
    kind          TEXT    NOT NULL,          -- lambda_reorg | uvvis_workflow | ...
    smiles        TEXT    NOT NULL,
    created_at    TEXT    NOT NULL,
    status        TEXT    NOT NULL,          -- pending | running | complete | error
    summary_json  TEXT                       -- derived quantities (λ_hole, λ_electron, Marcus rate, ...)
);

CREATE INDEX IF NOT EXISTS idx_workflows_status ON workflows(status);

-- Convenience view for the dashboard: jobs joined with their results.
CREATE VIEW IF NOT EXISTS v_jobs_with_results AS
SELECT
    j.id, j.signature, j.smiles, j.method, j.basis, j.charge, j.mult, j.job_kind,
    j.status, j.submitted_at, j.started_at, j.completed_at,
    j.log_path, j.ssh_jobid, j.workflow_id,
    r.homo_ev, r.lumo_ev, r.gap_ev, r.dipole_debye, r.energy_ev
FROM jobs j
LEFT JOIN results r ON r.job_id = j.id;
