"""Compose a correcting RM15 source with an encoded adaptive T instrument."""

from fractions import Fraction
from itertools import product
from math import comb

from ftqc_codes.gf2 import Echelon, dot
from ftqc_extractors.__main__ import digest, verify
from ftqc_extractors.algebra import bits, css_dimension, hgp, seed
from ftqc_extractors.magic import branches

from .distill import phase_sum
from .stochastic import probability, stochastic_result
from .supports import fields


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


def encoding(gate):
    n, hx, hz = hgp(*seed(gate["a"]), *seed(gate["b"]))
    k = css_dimension((n, hx, hz))
    x, z = gate["logical_x"], gate["logical_z"]
    if not dot(x, z) or any(dot(x, h) for h in hz) or any(dot(z, h) for h in hx):
        raise ValueError("invalid logical encoding port")
    ops = css_encoder(n, hx, x, z)
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


def instrument(native):
    # K[p,x] maps |b,a> to |b> for a=p XOR b, with amplitude
    # omega^(4xa+(2p+4x)b)/sqrt(2). Keep the classical transcript.
    basis = list(product(range(2), repeat=2))
    kraus = {}
    for p, m in basis:
        kraus[p, m] = {
            (b, a): (4 * m * a + (2 * p + 4 * m) * b) % 8
            for b, a in basis
            if a == p ^ b
        }
    for left, right in product(basis, repeat=2):
        terms = [
            row[right] - row[left]
            for row in kraus.values()
            if left[0] == right[0] and left in row and right in row
        ]
        if phase_sum(terms) != (2 * int(left == right), 0, 0, 0):
            raise ValueError("injection Kraus operators are not trace preserving")
    for resource_z in (0, 1):
        rows = native if not resource_z else branches(resource_z=1)
        if [(r["p"], r["x"]) for r in rows] != basis:
            raise ValueError("missing adaptive gate branch")
        for r in rows:
            p, m = r["p"], r["x"]
            cols = [
                [b, (phase + a + 4 * resource_z * a) % 8]
                for (b, a), phase in kraus[p, m].items()
            ]
            if (
                cols != r["columns"]
                or r["probability"] != "1/4"
                or (cols[1][1] - cols[0][1]) % 8 != 1 + 4 * resource_z
            ):
                raise ValueError(
                    "distilled T phase or feed-forward disagrees with consumer"
                )
    return {
        "transcripts": len(kraus),
        "completeness_matrix_entries": 16,
        "ideal_action": "T",
        "resource_Z_action": "Z T",
    }


def logical_t_result(spec):
    fields(
        spec,
        "kind",
        "factory",
        "handoff",
        "gate",
        "gate_receipt",
        "encoder_p",
        "online_p",
    )
    source = spec["factory"]
    if (
        source["kind"] != "local-stochastic"
        or source["source"]["kind"] != "level"
        or source["source"]["output_check"] != "distill15"
    ):
        raise ValueError("expected a checked 15-to-1 output distance")
    bound = stochastic_result(source)
    if bound["status"] != "proved" or not bound["also_bounds"]:
        raise ValueError("factory has an open output obligation")
    if spec["handoff"] != "fresh-product-css":
        raise ValueError("expected fresh resource and deterministic CSS encoding")
    gate, bank = spec["gate"], spec["gate_receipt"]
    if gate["kind"] != "t_gate" or gate["resource"] != "ideal-encoded-T":
        raise ValueError("expected an encoded logical T consumer")
    verify(gate, bank)
    result = bank["result"]
    if result["status"] != "proved":
        raise ValueError("consumer has no proved distance lower bound")
    adapter = encoding(gate)
    checked = instrument(result["ideal_branches"])
    rate, delta = probability(spec["online_p"]), Fraction(bound["upper"])
    encoder_p = probability(spec["encoder_p"])
    encoder_tail = min(Fraction(1), adapter["operation_count"] * encoder_p)
    input_delta = min(Fraction(1), delta + encoder_tail)
    cases = {}
    online = {}
    for name, case in result["cases"].items():
        n, d = case["locations"], case["evidence"]["lower"]
        online[name] = min(Fraction(1), comb(n, d) * rate**d)
        cases[name] = {
            "locations": n,
            "distance_lower": d,
            "online_fault_tail_upper": str(online[name]),
            "accepted_failure_upper": str(min(Fraction(1), input_delta + online[name])),
        }
    return {
        "status": "proved",
        "factory_output": "T|+>",
        "factory_p": bound["p"],
        "factory_noise_evidence": bound["noise_evidence"],
        "encoder": adapter,
        "instrument": checked,
        "resource_trace_distance_upper": str(delta),
        "fault_free_logical_T_half_diamond_upper": str(delta),
        "encoder_p": str(encoder_p),
        "encoder_any_fault_upper": str(encoder_tail),
        "encoded_resource_trace_distance_upper": str(input_delta),
        "online_p": str(rate),
        "cases": cases,
        "adaptive_readout_accepted_failure_upper": str(
            min(Fraction(1), input_delta + online["read_X"] + online["read_Y"])
        ),
        "scope": "fresh resource independent of data, reference and subsequent noise; ideal factory processing; clean encoded data; native detectors include ideal terminal syndromes",
        "encoding_noise_condition": "each listed preparation or CNOT has stochastic Pauli fault probability at most encoder_p; discard and frame processing ideal",
        "online_noise_condition": "stochastic Pauli faults with Pr[S contained in support] <= online_p^|S| on each potential arm before branch selection",
        "accepted_failure_metric": "joint probability of acceptance and failure of any test that passes with certainty on the correct ideal logical output; fixed-arm bounds are unnormalized",
    }
