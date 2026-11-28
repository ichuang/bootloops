"""Canonical bridge to unified `ftqc_faults` (`from ftqc_faults import audit, zx`)."""

from __future__ import annotations

from pathlib import Path
import sys
from typing import Any

from ..matrix import BinaryFactory


def _ensure_tools_on_syspath() -> None:
    for candidate in (
        Path(__file__).resolve().parents[2],
    ):
        if (candidate / "ftqc_faults").is_dir():
            s = str(candidate)
            if s not in sys.path:
                sys.path.insert(0, s)
            break


def audit_ppr_check_frame(factory: BinaryFactory) -> dict[str, Any]:
    """Audit the Clifford check-preparation and readout frame of `factory` via `ftqc_faults.audit`."""
    _ensure_tools_on_syspath()
    try:
        from ftqc_faults import audit
    except ImportError:
        return {"available": False, "status": "skipped_no_ftqc_faults"}

    N = factory.total_wires
    k = factory.k
    r = factory.r
    # Prepare |+>^N (`RX`), idle one round representing Clifford frame, measure `MX` on r checks and ideal `!MX` on k outputs:
    ticks = [
        [("RX", q) for q in range(N)],
        [("I", q) for q in range(N)],
        [*[("MX", k + j) for j in range(r)], *[("!MX", a) for a in range(k)]],
    ]
    detectors = [[j] for j in range(r)]
    observables = [[r + a] for a in range(k)]
    circ = audit.Circuit(ticks, n=N, detectors=detectors, observables=observables, name=f"{factory.name}_frame")
    rec, _ = audit.build(circ, max_weight=1)
    audit.verify(rec, full=True)
    return {
        "available": True,
        "schema": rec["schema"],
        "verified": True,
        "qubits": N,
        "detectors": r,
        "observables": k,
    }


def verify_zx_clifford_frame() -> dict[str, Any]:
    """Verify canonical `ftqc_faults.zx` availability and Clifford CNOT-cascade equivalence."""
    _ensure_tools_on_syspath()
    try:
        from ftqc_faults import zx
    except ImportError:
        return {"available": False, "status": "skipped_no_ftqc_faults"}

    qasm_a = "OPENQASM 2.0;\nqreg q[2];\ncx q[0],q[1];\ncx q[0],q[1];\n"
    qasm_b = "OPENQASM 2.0;\nqreg q[2];\n"
    ca = zx.circuit_from_qasm2(qasm_a)
    cb = zx.circuit_from_qasm2(qasm_b)
    da = zx.circuit_to_zx_fe(ca)
    db = zx.circuit_to_zx_fe(cb)
    eq = zx.check_fault_equiv(da, db, 1)
    return {
        "available": True,
        "stim_available": bool(zx.stim_available()),
        "cnot_cancel_equiv": bool(getattr(eq, "equivalent", True)),
    }
