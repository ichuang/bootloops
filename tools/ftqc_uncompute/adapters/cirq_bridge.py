"""Optional Cirq `Circuit` bridge (`skip-by-name` when `cirq` is absent)."""

from __future__ import annotations

from fractions import Fraction
import importlib.util
from typing import Any

from ..ir_kickmix import ANNOTATION_OPS, KickMixCircuit, KickMixOp


def cirq_available() -> bool:
    return importlib.util.find_spec("cirq") is not None


def to_cirq_circuit(circuit: KickMixCircuit) -> Any:
    """Convert an unconditioned or 1-level classically-controlled `KickMixCircuit` to `cirq.Circuit`."""
    if not cirq_available():
        raise ImportError("Optional dependency `cirq` is not installed")
    import cirq  # type: ignore[import-not-found]

    qs = cirq.LineQubit.range(circuit.num_qubits)
    cirq_ops: list[Any] = []
    for op in circuit.ops:
        if op.kind in ANNOTATION_OPS:
            continue
        k = op.kind
        gate_op: Any = None
        if k == "X":
            gate_op = cirq.X(qs[op.qubits[0]])
        elif k == "CX":
            gate_op = cirq.CNOT(qs[op.qubits[0]], qs[op.qubits[1]])
        elif k == "CCX":
            gate_op = cirq.TOFFOLI(qs[op.qubits[0]], qs[op.qubits[1]], qs[op.qubits[2]])
        elif k == "SWAP":
            gate_op = cirq.SWAP(qs[op.qubits[0]], qs[op.qubits[1]])
        elif k == "CSWAP":
            gate_op = cirq.CSWAP(qs[op.qubits[0]], qs[op.qubits[1]], qs[op.qubits[2]])
        elif k == "Z":
            gate_op = cirq.Z(qs[op.qubits[0]])
        elif k == "CZ":
            gate_op = cirq.CZ(qs[op.qubits[0]], qs[op.qubits[1]])
        elif k == "CCZ":
            gate_op = cirq.CCZ(qs[op.qubits[0]], qs[op.qubits[1]], qs[op.qubits[2]])
        elif k == "R":
            gate_op = cirq.reset(qs[op.qubits[0]])
        elif k == "HMR":
            q = qs[op.qubits[0]]
            key = f"b{op.bits[0]}"
            cirq_ops.append(cirq.H(q))
            cirq_ops.append(cirq.measure(q, key=key))
            cirq_ops.append(cirq.reset(q))
            continue
        if gate_op is not None:
            if op.condition is not None:
                gate_op = gate_op.with_classical_controls(f"b{op.condition}")
            cirq_ops.append(gate_op)
    return cirq.Circuit(cirq_ops)


def from_cirq_circuit(cirq_circ: Any) -> KickMixCircuit:
    if not cirq_available():
        raise ImportError("Optional dependency `cirq` is not installed")
    import cirq  # type: ignore[import-not-found]

    q_list = sorted(cirq_circ.all_qubits())
    q_map = {q: idx for idx, q in enumerate(q_list)}
    ops: list[KickMixOp] = []
    for cop in cirq_circ.all_operations():
        g = cop.gate
        qids = tuple(q_map[q] for q in cop.qubits)
        if g == cirq.X:
            ops.append(KickMixOp("X", qubits=qids))
        elif g == cirq.CNOT:
            ops.append(KickMixOp("CX", qubits=qids))
        elif g == cirq.TOFFOLI:
            ops.append(KickMixOp("CCX", qubits=qids))
        elif g == cirq.Z:
            ops.append(KickMixOp("Z", qubits=qids))
        elif g == cirq.CZ:
            ops.append(KickMixOp("CZ", qubits=qids))
        elif g == cirq.CCZ:
            ops.append(KickMixOp("CCZ", qubits=qids))
        elif g == cirq.S:
            ops.append(KickMixOp("Z_POW", qubits=qids, angle=Fraction(1, 2)))
        elif g == cirq.T:
            ops.append(KickMixOp("Z_POW", qubits=qids, angle=Fraction(1, 4)))
        else:
            raise ValueError(f"Unsupported Cirq operation {cop!r}")
    return KickMixCircuit(ops=tuple(ops), num_qubits=len(q_list))
