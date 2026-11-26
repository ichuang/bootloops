"""Physical outer coordinates for block-allocated gadgets."""

from .gadget import GadgetError, prep_layers


class ColumnError(GadgetError):
    def __init__(self, message, op=None):
        super().__init__(message)
        self.op = op


def layout(g):
    """Return wire coordinates and operations belonging to verified fresh specifications."""
    blocks = g.inputs + g.outputs
    if not blocks:
        raise ColumnError("physical columns require an interface block")
    n = len(blocks[0].qubits)
    refs = {q for b in g.inputs for q in b.refs}
    wires = [q for q in range(g.circuit.n) if q not in refs]
    if not n or len(wires) % n:
        raise ColumnError("qubits do not come in whole blocks")
    pos = {q: i % n for i, q in enumerate(wires)}
    if any([pos.get(q) for q in b.qubits] != list(range(n)) for b in blocks):
        raise ColumnError("interface order disagrees with physical columns")

    # Only a fresh encoded ancilla supplied as a specification can be prepared
    # across columns. Match its actual circuit, including the noisy boundary,
    # before exempting any operation. A declaration alone proves nothing.
    declared = [p for p in g.parts if p.get("method") == "spec"]
    exempt = set()
    for start in range(0, len(wires), n):
        if not declared:
            break
        qs = wires[start:start+n]
        touched = set(qs)
        first = next((k for k, tick in enumerate(g.circuit.ticks)
                      if any(touched.intersection(op.qs) for op in tick)), None)
        if first is None or first < g.body_start:
            continue
        for part in declared:
            code = next((b.code for b in blocks if b.code and b.code.name == part["code"]), None)
            if code is None or part["role"] not in ("zero", "plus"):
                continue
            layers = prep_layers(code, part["role"], qs)
            want = [[(False, *op) for op in layer] for layer in layers]
            want.append([(True, "I", q) for q in qs])
            ticks = g.circuit.ticks[first:first+len(want)]
            actual = [[op for op in t if touched.intersection(op.qs)] for t in ticks]
            if [sorted((op.noisy, op.name, *op.args()) for op in t) for t in actual] == list(map(sorted, want)):
                exempt.update(op for t in actual for op in t)
                declared.remove(part)
                break
    if declared:
        raise ColumnError("specified ancilla does not match a fresh preparation and its boundary")
    for tick in g.circuit.ticks[g.body_start:]:
        for op in tick:
            if op in exempt:
                continue
            qs = list(op.qs)
            if op.record is not None:
                qs.append(g.circuit.meas[op.record][2])
            if any(q not in pos for q in qs) or len({pos[q] for q in qs}) != 1:
                raise ColumnError(f"{op.text()} is not confined to a physical column (joins two positions or a reference)", op)
    return pos, exempt
