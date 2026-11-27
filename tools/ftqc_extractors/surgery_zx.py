"""Compile both ideal instruments of Theorem 9.1 and compare them with ftqc_zx."""

from types import SimpleNamespace
from ftqc_faults.zx.diagram import CircuitOp, CliffordCircuit
from ftqc_faults.zx.pauliweb import match_measurement_branches
from ftqc_faults.zx.equivalence import circuit_distance
from .algebra import bits, Limit
from .block_reading import memory, split_section, combine
from .surgery_sequence import append_cone
from .gross_probes import encoder
from .reading_audit import linear_frame


def analyze(s, budget):
    code = memory(s["memory"])
    n, hx, hz = code
    op = s["operation"]
    v, ax, az = memory(op["aux"])
    if (n + v + len(ax)) ** 2 > budget:
        raise Limit("ZX diagram width budget")
    _, measured = append_cone(code, op, n, list(hx), list(hz), s["max_degree"])
    f1 = op["f1"]
    inputs = tuple(range(n))

    def gate(name, *qubits):
        return CircuitOp(name, tuple(qubits))

    def circuit(width, gates, name):
        return CliffordCircuit(width, inputs, inputs, gates, name=name)

    gates = []
    for tick in encoder(SimpleNamespace(n=v, hx=ax, hz=az), "Z"):
        for kind, *sites in tick:
            gates.append(
                gate(
                    {"!RX": "PREP_PLUS", "!R": "PREP_0", "!CX": "CX"}[kind],
                    *(n + i for i in sites),
                )
            )
    for j, mask in enumerate(f1):
        gates.extend(gate("CX", q, n + j) for q in bits(mask))
    gates.extend(gate("MEAS_Z", n + j) for j in range(v))
    hom = circuit(n + v, gates, "homomorphic")
    edge = n + v
    gates = [gate("PREP_PLUS", edge + i) for i in range(len(ax))]
    for j, mask in enumerate(f1):
        gates.append(gate("PREP_0", n + j))
        gates.extend(gate("CX", q, n + j) for q in bits(mask))
        gates.extend(
            gate("CX", edge + i, n + j) for i, h in enumerate(ax) if h >> j & 1
        )
        gates.append(gate("MEAS_Z", n + j))
    gates.extend(gate("MEAS_X", edge + i) for i in range(len(ax)))
    surg = circuit(edge + len(ax), gates, "cone surgery")
    branch_map = match_measurement_branches(hom, surg)
    # Both vertex-record strings label the same coset modulo row(ax).
    # Fourier measurement of the edge workspace adds a phase Z(f1 * solution).
    frame = [combine(f1, col) for col in linear_frame(split_section(ax, v), len(ax))]
    return dict(
        status="proved",
        model="ideal Z measurement instrument up to split frame",
        measured=measured,
        vertex_coset_modulo=ax,
        split_Z_frame_columns=frame,
        zx_outcome_map=branch_map,
        physical_open_wire_checks={
            c.name: circuit_distance(c, w_max=1).to_dict() for c in (hom, surg)
        },
    )
