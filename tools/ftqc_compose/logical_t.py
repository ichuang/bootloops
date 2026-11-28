"""Compose a correcting RM15 source with an encoded adaptive T instrument."""

from fractions import Fraction
from itertools import product
from math import comb

from ftqc_gadgets.magic_transfer import css_encoder, verify_fresh
from ftqc_extractors.__main__ import verify
from ftqc_extractors.algebra import hgp, seed
from ftqc_extractors.magic import branches

from .distill import phase_sum
from .stochastic import probability, stochastic_result
from .supports import fields


def encoding(gate):
    n, hx, hz = hgp(*seed(gate["a"]), *seed(gate["b"]))
    x, z = gate["logical_x"], gate["logical_z"]
    return verify_fresh(n, hx, hz, x, z, css_encoder(n, hx, x, z))


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
    return _consume(spec, bound)


def _consume(spec, bound, *, locations=None):
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
    encoder_tail = min(Fraction(1), (adapter["operation_count"] if locations is None else locations) * encoder_p)
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
