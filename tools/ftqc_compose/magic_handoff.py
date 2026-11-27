"""Join a replayed source, a physical CSS transfer and the native T instrument."""

from fractions import Fraction

from ftqc_extractors.__main__ import digest
from ftqc_extractors.algebra import hgp, seed
from ftqc_gadgets.magic_transfer import fresh_receipt
from ftqc_magic.delivery import read_delivery

from .logical_t import _consume
from .supports import fields


def transfer_for(gate, port):
    n, hx, hz = hgp(*seed(gate["a"]), *seed(gate["b"]))
    if [gate["logical_x"], gate["logical_z"]] != port["pair"]:
        raise ValueError("destination logical pair differs from extractor")
    return fresh_receipt(n, hx, hz, port)


def consume_t(spec):
    fields(
        spec,
        "kind",
        "factory",
        "factory_receipt",
        "attempt",
        "gate",
        "gate_receipt",
        "port",
        "transfer_receipt",
        "encoder_p",
        "online_p",
        *(["encoded_origin"] if "encoded_origin" in spec else []),
    )
    attempt = spec["attempt"]
    if not isinstance(attempt, str) or not attempt.strip():
        raise ValueError("missing source attempt identity")
    source = read_delivery(spec["factory"], spec["factory_receipt"])
    if (
        source["state"] != "T"
        or source["groups"] != [[0]]
        or source["law"] != "exact-accepted-Z"
    ):
        raise ValueError("expected one verified T with exact Z law")
    transfer = transfer_for(spec["gate"], spec["port"])
    if transfer != spec["transfer_receipt"]:
        raise ValueError(
            "transfer receipt differs from the replayed port, circuit or frame"
        )
    delta = source["conditional_trace_distance_upper"]
    if "encoded_origin" in spec:
        origin = spec["encoded_origin"]
        if origin != dict(
            attempt=attempt,
            factory=source["producer"]["spec_sha256"],
            port=spec["port"],
            transfer=transfer,
            wire_order=list(range(transfer["input_wire"])),
        ):
            raise ValueError(
                "direct encoded origin differs from the replayed preparation"
            )
    # Both paths replay the original transfer once. Forwarding its encoded
    # output is an identity operation on the same ordered block.
    r = _consume(
        spec,
        dict(
            upper=delta,
            p=spec["factory_receipt"]["contract"]["input_error_rate"],
            noise_evidence=source["premises"],
        ),
        locations=len(transfer["operations"]),
    )
    acceptance = Fraction(source["acceptance_lower"])
    for case in r["cases"].values():
        case["attempt_joint_accepted_failure_upper"] = str(
            acceptance * Fraction(case["accepted_failure_upper"])
        )
    r.update(
        delivery=source,
        attempt=attempt,
        transfer=transfer,
        gate_receipt_sha256=digest(spec["gate_receipt"]),
        final_conditioning="joint acceptance and failure per attempted delivery; no final acceptance denominator inferred",
    )
    return r
