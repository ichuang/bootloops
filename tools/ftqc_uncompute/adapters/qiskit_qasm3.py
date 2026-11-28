"""Pure-stdlib OpenQASM 3.0 dynamic-circuit subset importer/exporter + optional Qiskit bridge.

Maps `HMR q[i] -> b[j]` to `h q[i]; b[j] = measure q[i]; reset q[i];` and classically
conditioned gates to `if (b[j]) { ... }`, using exact symbolic `pi * p / q` angles
(zero floating-point literals).
"""

from __future__ import annotations

from fractions import Fraction
import importlib.util
import re
from typing import Any

from ..ir_kickmix import ANNOTATION_OPS, Bit, KickMixCircuit, KickMixOp, Qubit, RegisterSpec


def qiskit_available() -> bool:
    return importlib.util.find_spec("qiskit") is not None


def _emit_qasm3_gate(op: KickMixOp) -> list[str]:
    k = op.kind
    qs = [f"q[{q}]" for q in op.qubits]
    bs = [f"b[{b}]" for b in op.bits]
    if k == "X":
        return [f"x {qs[0]};"]
    if k == "CX":
        return [f"cx {qs[0]}, {qs[1]};"]
    if k == "CCX":
        return [f"ccx {qs[0]}, {qs[1]}, {qs[2]};"]
    if k == "SWAP":
        return [f"swap {qs[0]}, {qs[1]};"]
    if k == "CSWAP":
        return [f"cswap {qs[0]}, {qs[1]}, {qs[2]};"]
    if k == "Z":
        return [f"z {qs[0]};"]
    if k == "CZ":
        return [f"cz {qs[0]}, {qs[1]};"]
    if k == "CCZ":
        return [f"h {qs[2]};", f"ccx {qs[0]}, {qs[1]}, {qs[2]};", f"h {qs[2]};"]
    if k == "Z_POW":
        assert op.angle is not None
        ang = op.angle
        if ang == Fraction(1, 2):
            return [f"s {qs[0]};"]
        if ang == Fraction(-1, 2) or ang == Fraction(3, 2):
            return [f"sdg {qs[0]};"]
        if ang == Fraction(1, 4):
            return [f"t {qs[0]};"]
        if ang == Fraction(-1, 4) or ang == Fraction(7, 4):
            return [f"tdg {qs[0]};"]
        return [f"p(pi * {ang.numerator} / {ang.denominator}) {qs[0]};"]
    if k == "NEG":
        return ["gphase(pi);"]
    if k == "R":
        return [f"reset {qs[0]};"]
    if k == "HMR":
        return [f"h {qs[0]};", f"{bs[0]} = measure {qs[0]};", f"reset {qs[0]};"]
    if k == "BIT_STORE0":
        return [f"{bs[0]} = 0;"]
    if k == "BIT_STORE1":
        return [f"{bs[0]} = 1;"]
    if k == "BIT_INVERT":
        return [f"{bs[0]} = !{bs[0]};"]
    return []


def kmx_to_qasm3(circuit: KickMixCircuit) -> str:
    """Export `circuit` to clean OpenQASM 3.0 with dynamic `if (b[k])` blocks."""
    lines: list[str] = [
        "OPENQASM 3.0;",
        'include "stdgates.inc";',
    ]
    for reg in circuit.registers:
        nm = reg.name if reg.name is not None else f"r{reg.reg_id}"
        w_list = ",".join(str(w.id) for w in reg.wires)
        lines.append(f"// @register r{reg.reg_id} {reg.kind} {nm} {w_list}")
    if circuit.num_qubits > 0:
        lines.append(f"qubit[{circuit.num_qubits}] q;")
    if circuit.num_bits > 0:
        lines.append(f"bit[{circuit.num_bits}] b;")

    indent = 0
    for op in circuit.ops:
        if op.kind in ANNOTATION_OPS:
            continue
        if op.kind == "PUSH_CONDITION":
            lines.append("  " * indent + f"if (b[{op.condition}]) {{")
            indent += 1
            continue
        if op.kind == "POP_CONDITION":
            indent = max(0, indent - 1)
            lines.append("  " * indent + "}")
            continue

        body = _emit_qasm3_gate(op)
        if op.condition is not None:
            lines.append("  " * indent + f"if (b[{op.condition}]) {{")
            for stmt in body:
                lines.append("  " * (indent + 1) + stmt)
            lines.append("  " * indent + "}")
        else:
            for stmt in body:
                lines.append("  " * indent + stmt)

    while indent > 0:
        indent -= 1
        lines.append("  " * indent + "}")
    return "\n".join(lines) + "\n"


def kmx_from_qasm3(text: str) -> KickMixCircuit:
    """Parse our dynamic-circuit OpenQASM 3.0 subset back into a `KickMixCircuit`,
    recognizing `h q[i]; b[j] = measure q[i]; reset q[i];` triples as native `HMR q[i] b[j]`.
    """
    nq = 0
    nb = 0
    ops: list[KickMixOp] = []
    registers: list[RegisterSpec] = []

    raw_stmts: list[str] = []
    for raw_line in text.splitlines():
        s_trim = raw_line.strip()
        if s_trim.startswith("// @register "):
            parts = s_trim[len("// @register ") :].split()
            if len(parts) == 4 and parts[0].startswith("r"):
                rid = int(parts[0][1:])
                rkind = parts[1]
                rname = parts[2]
                wids = [int(x) for x in parts[3].split(",") if x]
                wires = tuple(Qubit(w) if rkind == "q" else Bit(w) for w in wids)
                registers.append(RegisterSpec(reg_id=rid, name=rname, kind=rkind, wires=wires))
        line = raw_line.split("//", 1)[0].strip()
        if not line:
            continue
        raw_stmts.append(line)

    q_decl = re.compile(r"^qubit\[(\d+)\]\s+q;$")
    b_decl = re.compile(r"^bit\[(\d+)\]\s+b;$")
    if_hdr = re.compile(r"^if\s*\(\s*b\[(\d+)\]\s*\)\s*\{$")
    h_re = re.compile(r"^h\s+q\[(\d+)\];$")
    meas_re = re.compile(r"^b\[(\d+)\]\s*=\s*measure\s+q\[(\d+)\];$")
    rst_re = re.compile(r"^reset\s+q\[(\d+)\];$")
    q_idx_re = re.compile(r"q\[(\d+)\]")
    b_idx_re = re.compile(r"b\[(\d+)\]")
    p_re = re.compile(r"^p\(\s*pi\s*\*\s*(-?\d+)\s*/\s*(\d+)\s*\)\s+q\[(\d+)\];$")

    i = 0
    while i < len(raw_stmts):
        st = raw_stmts[i]
        if st.startswith(("OPENQASM", "include")):
            i += 1
            continue
        mq = q_decl.match(st)
        if mq:
            nq = max(nq, int(mq.group(1)))
            i += 1
            continue
        mb = b_decl.match(st)
        if mb:
            nb = max(nb, int(mb.group(1)))
            i += 1
            continue
        mif = if_hdr.match(st)
        if mif:
            ops.append(KickMixOp("PUSH_CONDITION", condition=int(mif.group(1))))
            i += 1
            continue
        if st == "}":
            ops.append(KickMixOp("POP_CONDITION"))
            i += 1
            continue

        # Recognize HMR triple: h q[x]; b[y] = measure q[x]; reset q[x];
        if i + 2 < len(raw_stmts):
            mh = h_re.match(st)
            mm = meas_re.match(raw_stmts[i + 1])
            mr = rst_re.match(raw_stmts[i + 2])
            if mh and mm and mr and mh.group(1) == mm.group(2) == mr.group(1):
                ops.append(KickMixOp("HMR", qubits=(int(mh.group(1)),), bits=(int(mm.group(1)),)))
                i += 3
                continue
            # Recognize CCZ triple: h q[z]; ccx q[x], q[y], q[z]; h q[z];
            mh2 = h_re.match(raw_stmts[i + 2])
            if mh and mh2 and mh.group(1) == mh2.group(1) and raw_stmts[i + 1].startswith("ccx "):
                qs_ccz = [int(x) for x in q_idx_re.findall(raw_stmts[i + 1])]
                if len(qs_ccz) == 3 and qs_ccz[2] == int(mh.group(1)):
                    ops.append(KickMixOp("CCZ", qubits=tuple(qs_ccz)))
                    i += 3
                    continue

        mp = p_re.match(st)
        if mp:
            ang = Fraction(int(mp.group(1)), int(mp.group(2)))
            ops.append(KickMixOp("Z_POW", qubits=(int(mp.group(3)),), angle=ang))
            i += 1
            continue

        if st == "gphase(pi);":
            ops.append(KickMixOp("NEG"))
        elif st.endswith("= 0;"):
            bs = [int(x) for x in b_idx_re.findall(st)]
            ops.append(KickMixOp("BIT_STORE0", bits=(bs[0],)))
        elif st.endswith("= 1;"):
            bs = [int(x) for x in b_idx_re.findall(st)]
            ops.append(KickMixOp("BIT_STORE1", bits=(bs[0],)))
        elif "= !" in st:
            bs = [int(x) for x in b_idx_re.findall(st)]
            ops.append(KickMixOp("BIT_INVERT", bits=(bs[0],)))
        else:
            head = st.rstrip(";").split(None, 1)[0].lower()
            qs = tuple(int(x) for x in q_idx_re.findall(st))
            simple_map = {
                "x": ("X", None),
                "cx": ("CX", None),
                "ccx": ("CCX", None),
                "swap": ("SWAP", None),
                "cswap": ("CSWAP", None),
                "z": ("Z", None),
                "cz": ("CZ", None),
                "reset": ("R", None),
                "s": ("Z_POW", Fraction(1, 2)),
                "sdg": ("Z_POW", Fraction(-1, 2)),
                "t": ("Z_POW", Fraction(1, 4)),
                "tdg": ("Z_POW", Fraction(-1, 4)),
            }
            if head not in simple_map:
                raise ValueError(f"Unsupported OpenQASM 3 statement: {st!r}")
            kind, ang = simple_map[head]
            ops.append(KickMixOp(kind, qubits=qs, angle=ang))
        i += 1

    return KickMixCircuit(ops=tuple(ops), registers=tuple(registers), num_qubits=nq, num_bits=nb)


def to_qiskit_circuit(circuit: KickMixCircuit) -> Any:
    """Convert `circuit` to `qiskit.QuantumCircuit` via OpenQASM 3 if `qiskit.qasm3` is installed."""
    if not qiskit_available():
        raise ImportError("Optional dependency `qiskit` is not installed")
    import qiskit.qasm3  # type: ignore[import-not-found]

    return qiskit.qasm3.loads(kmx_to_qasm3(circuit))
