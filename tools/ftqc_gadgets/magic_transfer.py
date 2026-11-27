"""CSS isometries for delivery of an open resource port."""

from ftqc_codes.gf2 import Echelon, dot
from ftqc_extractors.algebra import bits, css_dimension
from ftqc_extractors.__main__ import digest


def css_encoder(n, hx, logical_x, logical_z):
    """Prepare the positive X-check superposition, copy in bX, erase input b."""
    group = Echelon(hx)
    pivots = set(group.pivots)
    ops = [("RX" if q in pivots else "R", q) for q in range(n)]
    for p, row in zip(group.pivots, group.basis, strict=True):
        ops.extend(("CX", p, q) for q in bits(row ^ (1 << p)))
    ops.extend(("CX", n, q) for q in bits(logical_x))
    ops.extend(("CX", q, n) for q in bits(logical_z))
    return ops


def verify_fresh(n, hx, hz, x, z, ops):
    if (
        type(n) is not int
        or n < 1
        or any(type(v) is not int or not 0 <= v < 1 << n for v in [*hx, *hz])
        or any(dot(a, b) for a in hx for b in hz)
    ):
        raise ValueError("invalid CSS checks")
    k = css_dimension((n, hx, hz))
    if not dot(x, z) or any(dot(x, h) for h in hz) or any(dot(z, h) for h in hx):
        raise ValueError("invalid logical encoding port")
    # Each wire carries a binary linear form in b (bit zero) and the
    # positive-superposition variables from RX. CX introduces no phases.
    if (
        [op[1] for op in ops[:n]] != list(range(n))
        or any(op[0] not in ("R", "RX") for op in ops[:n])
        or any(op[0] != "CX" for op in ops[n:])
    ):
        raise ValueError("encoder must initialize the memory once, followed by CNOTs")
    forms = [0] * n + [1]
    variables = 0
    for op, *qs in ops:
        if op == "RX":
            variables += 1
            forms[qs[0]] = 1 << variables
        elif op == "CX":
            c, t = qs
            forms[t] ^= forms[c]
        elif op != "R":
            raise ValueError("unsupported CSS encoder operation")
    columns = [
        sum(((v >> j) & 1) << q for q, v in enumerate(forms[:n]))
        for j in range(variables + 1)
    ]
    span = Echelon(hx)
    if (
        forms[n]
        or columns[0] != x
        or variables != span.rank
        or Echelon(columns[1:]).rank != variables
        or any(c not in span for c in columns[1:])
    ):
        raise ValueError(
            "encoder does not realize the two logical cosets and erase its input"
        )
    # A Z fault on the raw input after its first fanout passes
    # through the remaining erasure CNOTs as the full encoded logical Z.
    site = next(i for i, op in enumerate(ops) if op[0] == "CX" and op[1] == n)
    error = 1 << n
    for _, control, target in ops[site + 1 :]:
        if error >> target & 1:
            error ^= 1 << control
    if error != (1 << n) | z:
        raise ValueError("encoder logical Z witness failed forward transport")
    return {
        "raw_to_gate_accepted_fault_distance": 1,
        "logical_Z_witness": {
            "location": site,
            "op": list(ops[site]),
            "after": "ZI",
            "output": "logical Z on resource; Z on discarded zero wire",
        },
        "memory_qubits": n,
        "logical_qubits": k,
        "logical_x": x,
        "logical_z": z,
        "spectators": k - 1,
        "positive_superposition_rank": variables,
        "operation_count": len(ops),
        "operations_sha256": digest(ops),
        "ideal_map": "deterministic CSS encoding, resource spectators initialized to logical zero in a Z basis commuting with logical X",
    }


def fresh_receipt(n, hx, hz, port):
    if (
        set(port) != {"block", "pair", "frame", "correction"}
        or not isinstance(port["block"], str)
        or not port["block"]
    ):
        raise ValueError("expected a named fresh block and signed logical pair")
    x, z = port["pair"]
    if any(type(v) is not int or not 0 < v < 1 << n for v in (x, z)):
        raise ValueError("logical representative outside destination")
    frame, correction = port["frame"], port["correction"]
    if frame not in ("I", "Z") or correction != frame:
        raise ValueError("unproved frame correction")
    ops = css_encoder(n, hx, x, z)
    checked = verify_fresh(n, hx, hz, x, z, ops)
    # Coherent image of (|0,R0> + |1,R1>)/sqrt(2): the verifier
    # checks both disjoint CSS cosets, identical positive amplitudes and
    # erasure. This fixes the off-diagonal blocks for arbitrary references.
    frames = [["Z", n], ["Z", n]] if frame == "Z" else []
    return dict(
        schema="ftqc-gadgets-fresh-transfer-1",
        port=port,
        checks=[list(hx), list(hz)],
        input_wire=n,
        discarded=[n],
        frame_ops=frames,
        operations=frames + [list(op) for op in ops],
        bell_isometry="positive equal CSS cosets correlated coherently with the untouched reference",
        verification=checked,
    )


def transfer_receipt(spec):
    if spec.get("kind") == "fresh" and set(spec) == {"kind", "n", "hx", "hz", "port"}:
        return fresh_receipt(spec["n"], spec["hx"], spec["hz"], spec["port"])
    if spec.get("kind") == "swap" and set(spec) == {"kind", "block", "slot"}:
        from .occupied_transfer import occupied_swap

        return occupied_swap(spec["block"], spec["slot"])
    raise ValueError("expected a fresh or swap transfer specification")


def verify_transfer(rec):
    if rec["schema"] == "ftqc-gadgets-fresh-transfer-1":
        hx, hz = rec["checks"]
        expected = fresh_receipt(rec["input_wire"], hx, hz, rec["port"])
    elif rec["schema"] == "ftqc-gadgets-occupied-swap-1":
        from .occupied_transfer import occupied_swap

        expected = occupied_swap(
            dict(n=rec["raw"], checks=rec["checks"], basis=rec["basis"]), rec["slot"]
        )
    else:
        raise ValueError("unknown transfer receipt")
    if digest(rec) != digest(expected):
        raise ValueError("transfer differs from replay")
