"""Cache lookup by signature hash.

signature = sha1(canonical_smiles + method + basis + charge + mult).hexdigest()[:16]

Cache hit = a `jobs` row exists with the same signature and status='complete'.
The function returns the result row so callers can short-circuit g16.
"""

from __future__ import annotations

import hashlib


def signature(
    canonical_smiles: str,
    method: str,
    basis: str,
    charge: int,
    mult: int,
    *,
    job_kind: str = "properties",
    extras: str = "",
) -> str:
    payload = (
        f"{canonical_smiles}|{method.lower()}|{basis.lower()}|{charge}|{mult}|{job_kind}|{extras}"
    )
    return hashlib.sha1(payload.encode()).hexdigest()[:16]
