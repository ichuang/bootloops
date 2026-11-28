"""Bind an inner gadget's signed logical instrument to an outer operation."""

from ftqc_faults.audit.circuit import Circuit
from ftqc_faults.audit.stab import Sim, _g_sum
from ftqc_codes import gf2

from .certify import Refusal
from .gadget import logical_x, logical_z

ARITY = {"cnot": (2, 2), "idle": (1, 1), "ec": (1, 1), "h": (1, 1), "s": (1, 1),
         "measure_z": (1, 0), "measure_x": (1, 0), "prep_zero": (0, 1), "prep_plus": (0, 1)}


def physical(mask, block):
    return sum(1 << block.qubits[i] for i in gf2.supp(mask))


def logical_action(g, role, code):
    if (len(g.inputs), len(g.outputs)) != ARITY[role]:
        raise Refusal(f"{role}: wrong logical input/output arity")
    blocks = g.inputs + g.outputs
    if any(not b.code or (b.code.n, b.code.hx, b.code.hz) != (code.n, code.hx, code.hz)
           or len(b.qubits) != code.n or len(set(b.qubits)) != code.n for b in blocks):
        raise Refusal(f"{role}: every port must use the same code in its declared wire order")
    if any(len(b.refs) != 1 for b in g.inputs):
        raise Refusal(f"{role}: each input needs a single logical reference")
    lx, = logical_x(code)
    lz, = logical_z(code)
    if not gf2.dot(lx, lz):
        raise Refusal(f"{role}: logical X and Z do not anticommute")
    refs = [1 << b.refs[0] for b in g.inputs]
    ref_wires = {b.refs[0] for b in g.inputs}
    for ports in (g.inputs, g.outputs):
        wires = [q for b in ports for q in b.qubits]
        if len(set(wires)) != len(wires) or ref_wires.intersection(wires) or len(ref_wires) != len(refs):
            raise Refusal(f"{role}: overlapping logical interfaces")
    if any(ref_wires.intersection(op.qs) for tick in g.circuit.ticks[g.body_start:] for op in tick):
        raise Refusal(f"{role}: body acts on an input reference")

    pre = Sim(g.circuit.n).run(Circuit(g.circuit.ticks[:g.body_start], g.circuit.n))
    sim = Sim(g.circuit.n).run(g.circuit)
    frames = {}

    def sign(state, x, z, label):
        try:
            return state._express(x, z)[0]
        except AssertionError:
            raise Refusal(f"{role}: logical relation is absent for {label}") from None

    def relation(x, z, label):
        value = sign(sim, x, z, label)
        if value & 1:
            raise Refusal(f"{role}: wrong sign for {label}")
        records = gf2.supp(sim.frame_mask(value))
        if records:
            frames[label] = records
        return records

    for r, block in zip(refs, g.inputs):
        if (sign(pre, r ^ physical(lx, block), 0, "input X")
                or sign(pre, 0, r ^ physical(lz, block), "input Z")):
            raise Refusal(f"{role}: input is not the canonical encoded Bell pair")
    for i, block in enumerate(g.outputs):
        for j, row in enumerate(code.hx):
            relation(physical(row, block), 0, f"block {i} X check {j}")
        for j, row in enumerate(code.hz):
            relation(0, physical(row, block), f"block {i} Z check {j}")
    xs = [physical(lx, b) for b in g.outputs]
    zs = [physical(lz, b) for b in g.outputs]
    if role in ("measure_x", "measure_z"):
        records = relation(refs[0] if role == "measure_x" else 0,
                           refs[0] if role == "measure_z" else 0, "measurement")
        if not records:
            raise Refusal(f"{role}: logical measurement outcome is missing")
    elif role == "prep_zero":
        relation(0, zs[0], "prepared Z")
    elif role == "prep_plus":
        relation(xs[0], 0, "prepared X")
    else:
        # Choi generators for P on a reference and U P U† on the outputs.
        if role == "cnot":
            xs[0] ^= xs[1]
            zs[1] ^= zs[0]
        for i, r in enumerate(refs):
            if role == "h":
                relation(r, zs[i], "input X")
                relation(xs[i], r, "input Z")
            elif role == "s":
                want_y_sign = ((_g_sum(lx, 0, 0, lz) + 1) // 2) & 1
                val = sign(sim, r ^ xs[i], zs[i], "input X")
                if (val & 1) != want_y_sign:
                    raise Refusal(f"{role}: wrong sign for input X")
                recs = gf2.supp(sim.frame_mask(val))
                if recs:
                    frames["input X"] = recs
                relation(0, r ^ zs[i], "input Z")
            else:
                relation(r ^ xs[i], 0, f"input {i} X")
                relation(0, r ^ zs[i], f"input {i} Z")
    return {"operation": role, "record_parities": frames}
