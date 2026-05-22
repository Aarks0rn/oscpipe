"""LocalBackend unit tests — no g16 required.

submit() invokes a real g16 binary and is covered by test_real_workstation.py.
All other behaviour (poll, fetch_log, preflight) is tested in isolation.
"""

from pathlib import Path

from oscpipe.dispatch.local import LocalBackend

_NORMAL = " Normal termination of Gaussian 16 at ...\n"
_ERROR = " Error termination via Lnk1e in /...\n"


def _backend(tmp_path):
    return LocalBackend(str(tmp_path / "work"))


def _write_log(backend, label, content):
    log = backend.work_dir / f"{label}.log"
    log.write_text(content)


def test_localbackend_poll_unknown(tmp_path):
    assert _backend(tmp_path).poll("nope") == "unknown"


def test_localbackend_poll_complete(tmp_path):
    b = _backend(tmp_path)
    _write_log(b, "h2", _NORMAL)
    assert b.poll("h2") == "complete"


def test_localbackend_poll_error(tmp_path):
    b = _backend(tmp_path)
    _write_log(b, "h2", _ERROR)
    assert b.poll("h2") == "error"


def test_localbackend_fetch_log_copies_to_outdir(tmp_path):
    b = _backend(tmp_path)
    _write_log(b, "h2", _NORMAL)
    out_dir = tmp_path / "out"
    path = b.fetch_log("h2", "h2", str(out_dir))
    assert path == str(out_dir / "h2.log")
    assert "Normal termination" in Path(path).read_text()


def test_localbackend_fetch_log_same_dir_no_copy(tmp_path):
    b = _backend(tmp_path)
    _write_log(b, "h2", _NORMAL)
    path = b.fetch_log("h2", "h2", str(b.work_dir))
    assert path == str(b.work_dir / "h2.log")


def test_localbackend_preflight_no_g16(tmp_path):
    b = LocalBackend(str(tmp_path / "work"), exe="g16-does-not-exist")
    checks = {name: ok for name, ok, _ in b.preflight()}
    assert checks["g16"] is False
    assert checks["work_dir"] is True


def test_localbackend_preflight_work_dir_writable(tmp_path):
    b = _backend(tmp_path)
    checks = {name: ok for name, ok, _ in b.preflight()}
    assert checks["work_dir"] is True
