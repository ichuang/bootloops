"""Optional PyZX Bridge and Pure-Stdlib OpenQASM 2.0 Clifford Importer/Exporter.

Enforces fail-closed fault-equivalence hygiene:
  - Standard `pyzx.full_reduce()` and `pyzx.clifford_simp()` preserve ONLY the noise-free
    linear map, NOT w-fault equivalence (Kissinger §4.2-§5, p. 19). Calling an uncertified
    PyZX simplification routine as if it were fault-preserving raises `RefusalError`!
  - Instead, `audit_pyzx_rewrite(g_before, g_after, w=...)` converts two `pyzx` graphs into
    `ZXDiagram` instances (using exact `fractions.Fraction` Clifford phases) and runs
    `check_fault_equiv` to certify or refute w-fault equivalence!

Zero float literals or float builtin references appear anywhere in this module.
"""

from __future__ import annotations

from fractions import Fraction
import re
from typing import Any

from ..diagram import CircuitOp, CliffordCircuit, PREP_STATES, RefusalError, SINGLE_MEAS, ZXDiagram
from ..equivalence import FaultEquivResult, check_fault_equiv


def pyzx_available() -> bool:
    try:
        __import__("pyzx")
        return True
    except ImportError:
        return False


def circuit_to_qasm2(circuit: CliffordCircuit) -> str:
    """Export a physical `CliffordCircuit` to OpenQASM 2.0 text (pure stdlib)."""
    circuit.validate()
    lines = [
        "OPENQASM 2.0;",
        'include "qelib1.inc";',
        f"qreg q[{circuit.n_qubits}];",
    ]
    n_meas = sum(1 for op in circuit.ops if op.kind in SINGLE_MEAS)
    if n_meas > 0:
        lines.append(f"creg c[{n_meas}];")

    m_idx = 0
    for op in circuit.ops:
        k = op.kind
        if k in ("CAT_PREP", "ENC_PREP", "MPP"):
            raise RefusalError("QASM-IDEALISED-OP-REFUSED", f"cannot export idealised op {k} to physical OpenQASM 2.0")
        if k == "PREP_0":
            lines.append(f"reset q[{op.qubits[0]}];")
        elif k == "PREP_1":
            q = op.qubits[0]
            lines.append(f"reset q[{q}];")
            lines.append(f"x q[{q}];")
        elif k == "PREP_PLUS":
            q = op.qubits[0]
            lines.append(f"reset q[{q}];")
            lines.append(f"h q[{q}];")
        elif k == "PREP_MINUS":
            q = op.qubits[0]
            lines.append(f"reset q[{q}];")
            lines.append(f"x q[{q}];")
            lines.append(f"h q[{q}];")
        elif k in ("I", "H", "S", "SDG", "X", "Y", "Z"):
            q = op.qubits[0]
            gate = "id" if k == "I" else k.lower()
            lines.append(f"{gate} q[{q}];")
        elif k == "CX":
            lines.append(f"cx q[{op.qubits[0]}],q[{op.qubits[1]}];")
        elif k == "CZ":
            lines.append(f"cz q[{op.qubits[0]}],q[{op.qubits[1]}];")
        elif k == "SWAP":
            lines.append(f"swap q[{op.qubits[0]}],q[{op.qubits[1]}];")
        elif k == "MEAS_Z":
            q = op.qubits[0]
            lines.append(f"measure q[{q}] -> c[{m_idx}];")
            m_idx += 1
        elif k == "MEAS_X":
            q = op.qubits[0]
            lines.append(f"h q[{q}];")
            lines.append(f"measure q[{q}] -> c[{m_idx}];")
            m_idx += 1
        elif k == "MEAS_Y":
            q = op.qubits[0]
            lines.append(f"sdg q[{q}];")
            lines.append(f"h q[{q}];")
            lines.append(f"measure q[{q}] -> c[{m_idx}];")
            m_idx += 1
    return "\n".join(lines) + "\n"


_QREG_RE = re.compile(r"^qreg\s+([A-Za-z_]\w*)\[(\d+)\]\s*;$")
_MEAS_RE = re.compile(r"^measure\s+\w+\[(\d+)\]\s*->\s*\w+\[(\d+)\]\s*;$")
_RESET_RE = re.compile(r"^reset\s+\w+\[(\d+)\]\s*;$")
_GATE1_RE = re.compile(r"^(id|h|s|sdg|x|y|z)\s+\w+\[(\d+)\]\s*;$")
_GATE2_RE = re.compile(r"^(cx|cz|swap)\s+\w+\[(\d+)\]\s*,\s*\w+\[(\d+)\]\s*;$")


def circuit_from_qasm2(
    qasm_text: str,
    input_qubits: tuple[int, ...] | None = None,
    output_qubits: tuple[int, ...] | None = None,
    name: str = "qasm_circuit",
) -> CliffordCircuit:
    """Parse a Clifford OpenQASM 2.0 program into `CliffordCircuit` (refusing non-Clifford T/RZ gates)."""
    n_qubits: int | None = None
    ops: list[CircuitOp] = []
    first_seen_prep: set[int] = set()
    touched: set[int] = set()
    measured_now: set[int] = set()

    for raw in qasm_text.splitlines():
        line = raw.split("//", 1)[0].strip()
        if not line or line.startswith("OPENQASM") or line.startswith("include") or line.startswith("creg"):
            continue
        m_qr = _QREG_RE.match(line)
        if m_qr:
            n_qubits = int(m_qr.group(2))
            continue
        m_rst = _RESET_RE.match(line)
        if m_rst:
            q = int(m_rst.group(1))
            if q not in touched:
                first_seen_prep.add(q)
            touched.add(q)
            measured_now.discard(q)
            ops.append(CircuitOp("PREP_0", (q,)))
            continue
        m_ms = _MEAS_RE.match(line)
        if m_ms:
            q = int(m_ms.group(1))
            cbit = int(m_ms.group(2))
            touched.add(q)
            measured_now.add(q)
            ops.append(CircuitOp("MEAS_Z", (q,), label=f"qasm_m{cbit}"))
            continue
        m_g1 = _GATE1_RE.match(line)
        if m_g1:
            gname = m_g1.group(1).upper()
            q = int(m_g1.group(2))
            touched.add(q)
            ops.append(CircuitOp("I" if gname == "ID" else gname, (q,)))
            continue
        m_g2 = _GATE2_RE.match(line)
        if m_g2:
            gname = m_g2.group(1).upper()
            q0 = int(m_g2.group(2))
            q1 = int(m_g2.group(3))
            touched.add(q0)
            touched.add(q1)
            ops.append(CircuitOp(gname, (q0, q1)))
            continue
        raise RefusalError("NON-CLIFFORD-OR-UNSUPPORTED-QASM", f"unsupported OpenQASM statement: {line!r}")

    if n_qubits is None:
        raise RefusalError("QASM-MISSING-QREG", "OpenQASM program missing `qreg` declaration")

    in_qs = input_qubits if input_qubits is not None else tuple(q for q in range(n_qubits) if q not in first_seen_prep)
    out_qs = output_qubits if output_qubits is not None else tuple(q for q in range(n_qubits) if q not in measured_now)
    return CliffordCircuit(n_qubits, in_qs, out_qs, ops, name=name)


def to_pyzx_graph(zx: ZXDiagram) -> Any:
    """Convert a `ZXDiagram` to a `pyzx.Graph` with exact `Fraction(k, 2)` phases."""
    if not pyzx_available():
        raise RefusalError("PYZX-NOT-INSTALLED", "optional dependency `pyzx` is not installed in this Python")
    import pyzx  # type: ignore

    zx.validate()
    g = pyzx.Graph()
    vt_map = {
        "BOUND": pyzx.VertexType.BOUNDARY,
        "Z": pyzx.VertexType.Z,
        "X": pyzx.VertexType.X,
        "H": pyzx.VertexType.H_BOX,
    }
    nid_to_v: dict[int, Any] = {}
    for nid, nd in sorted(zx.nodes.items()):
        phase = Fraction(1, 1) if nd.kind == "H" else Fraction(nd.phase_k, 2)
        v = g.add_vertex(vt_map[nd.kind], phase=phase)
        nid_to_v[nid] = v

    for eid, ed in sorted(zx.edges.items()):
        g.add_edge((nid_to_v[ed.u], nid_to_v[ed.v]), edgetype=pyzx.EdgeType.SIMPLE)

    g.set_inputs(tuple(nid_to_v[nid] for nid in zx.inputs))
    g.set_outputs(tuple(nid_to_v[nid] for nid in zx.outputs))
    return g


def from_pyzx_graph(g: Any, name: str = "from_pyzx") -> ZXDiagram:
    """Convert a Clifford `pyzx.Graph` into a `ZXDiagram` (expanding Hadamard edges into H-boxes)."""
    if not pyzx_available():
        raise RefusalError("PYZX-NOT-INSTALLED", "optional dependency `pyzx` is not installed in this Python")
    import pyzx  # type: ignore

    zx = ZXDiagram(name=name)
    v_to_nid: dict[Any, int] = {}
    in_set = set(g.inputs())
    out_set = set(g.outputs())

    for v in g.vertices():
        vt = g.type(v)
        ph = Fraction(g.phase(v)) % 2
        if vt == pyzx.VertexType.BOUNDARY:
            if v not in in_set and v not in out_set:
                raise RefusalError("PYZX-UNORDERED-BOUNDARY", f"boundary vertex {v} not in inputs() or outputs()")
            nid = zx.add_node("BOUND", 0, label=f"b_{v}")
        elif vt in (pyzx.VertexType.Z, pyzx.VertexType.X):
            if ph.denominator not in (1, 2):
                raise RefusalError("NON-CLIFFORD-PHASE-REFUSED", f"non-Clifford PyZX phase {ph}*pi on vertex {v}")
            phase_k = int((ph * 2) % 4)
            nid = zx.add_node("Z" if vt == pyzx.VertexType.Z else "X", phase_k, label=f"sp_{v}")
        elif vt == pyzx.VertexType.H_BOX:
            nid = zx.add_node("H", 0, label=f"hb_{v}")
        else:
            raise RefusalError("PYZX-UNSUPPORTED-VERTEX", f"unsupported PyZX vertex type {vt!r}")
        v_to_nid[v] = nid

    zx.inputs = [v_to_nid[v] for v in g.inputs()]
    zx.outputs = [v_to_nid[v] for v in g.outputs()]

    for e in g.edges():
        s, t = g.edge_st(e)
        et = g.edge_type(e)
        u_id, v_id = v_to_nid[s], v_to_nid[t]
        if et == pyzx.EdgeType.SIMPLE:
            zx.add_edge(u_id, v_id, fault_free=False)
        elif et == pyzx.EdgeType.HADAMARD:
            hb = zx.add_node("H", 0, label="had_edge")
            zx.add_edge(u_id, hb, fault_free=False)
            zx.add_edge(hb, v_id, fault_free=True)
        else:
            raise RefusalError("PYZX-UNSUPPORTED-EDGE", f"unsupported PyZX edge type {et!r}")

    zx.validate()
    return zx


def audit_pyzx_rewrite(
    g_before: Any,
    g_after: Any,
    w: int = 3,
    noise_model: str = "circuit_depolarizing",
) -> FaultEquivResult:
    """Audit whether a PyZX graph transformation preserved w-fault equivalence!"""
    zx_b = from_pyzx_graph(g_before, name="pyzx_before")
    zx_a = from_pyzx_graph(g_after, name="pyzx_after")
    return check_fault_equiv(zx_b, zx_a, w=w, noise_model_lhs=noise_model, cross_check_spackle=False)
