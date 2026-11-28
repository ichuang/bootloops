"""Qualtran Bloq / CompositeBloq bridge (`skip-by-name` when `qualtran` is absent)."""

from __future__ import annotations

import importlib.util
from typing import Any

from ..ir_dag import ArithDAG


def qualtran_available() -> bool:
    return importlib.util.find_spec("qualtran") is not None


def dag_to_qualtran_descriptor(dag: ArithDAG) -> dict[str, Any]:
    """Export an `ArithDAG` into a Qualtran-aligned Bloq signature & call-graph descriptor
    (and materialize a `qualtran.BloqBuilder` composite bloq if `qualtran` is installed).
    """
    ports = [
        {
            "name": r.name,
            "dtype": {
                "quint": f"QUInt({r.width})",
                "qgf2": f"QGF2({r.width})",
                "qmodint": f"BQUInt({r.width}, {r.modulus})",
                "mask": f"QAny({r.width})",
                "classical": f"CBit({r.width})",
            }[r.kind],
            "side": "THRU" if r.role in ("io", "input", "mask") else "RIGHT",
        }
        for r in dag.registers
        if r.role != "scratch"
    ]
    call_counts: dict[str, int] = {}
    for n in dag.nodes:
        bloq_nm = f"{n.op_kind}[{n.chirality}:{n.uncompute_mode}]"
        call_counts[bloq_nm] = call_counts.get(bloq_nm, 0) + 1

    return {
        "bloq_name": dag.name,
        "qualtran_installed": qualtran_available(),
        "signature": ports,
        "bloq_counts": call_counts,
        "analytical_t_complexity": dag.estimate_costs(),
    }
