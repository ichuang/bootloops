"""Exact accepted batch distributions and an ideal logical S consumer."""

from fractions import Fraction


def s_injection():
    # Data controls a CNOT onto (|0> + sign*i|1>)/sqrt(2), followed by
    # Z readout m of the resource and the correction Z**m on data.
    for wrong in (0, 1):
        phase = 1j * (-1) ** wrong
        for m in (0, 1):
            column = [(phase if d ^ m else 1) * (-1) ** (m * d) for d in (0, 1)]
            if column[1] != column[0] * phase:
                raise ValueError("incorrect logical S injection branch")
    return "ideal decoded Y resource, CNOT, Z readout and outcome-controlled Z; wrong Y sign gives Z after S"


def result(s):
    from ftqc_extractors.injection_export import y_law
    from ftqc_extractors.algebra import Limit

    if set(s) != {"kind", "case", "locations", "frame", "slots", "state_budget"}:
        raise ValueError("unsupported batch contract fields")
    try:
        spec, evidence, accepted, channel_count = y_law(s)
    except Limit as exc:
        raise ValueError(str(exc)) from exc
    k = len(accepted).bit_length() - 1
    a = sum(accepted)
    if not a:
        raise ValueError("zero acceptance")
    distribution = [p / a for p in accepted]
    slots = s["slots"]
    if (
        not isinstance(slots, list)
        or not slots
        or any(type(i) is not int or not 0 <= i < k for i in slots)
        or len(set(slots)) != len(slots)
    ):
        raise ValueError("the consumer requires distinct output slots")
    consumer = [Fraction(0)] * (1 << len(slots))
    for m, p in enumerate(distribution):
        index = sum(((m >> slot) & 1) << j for j, slot in enumerate(slots))
        consumer[index] += p
    return dict(
        status="proved",
        kind="accepted simultaneous Y source and ideal S consumers",
        scope="independent listed physical channels; all unlisted locations ideal; ideal terminal decoding and transfer",
        noise_scope="all physical channels"
        if s["locations"] == "all"
        else "selected physical channels",
        circuit_sha256=evidence["circuit_sha256"],
        source_sha256=evidence["source_sha256"],
        physical_channels=channel_count,
        raw_sites=spec["preparation"]["carriers"],
        frame=s["frame"],
        acceptance=str(a),
        joint_error_masks=[str(p) for p in distribution],
        accepted_states_per_attempt=str(k * a),
        qubit_layers_per_accepted_state=str(
            evidence["cost"]["physical_qubits"]
            * evidence["cost"]["two_qubit_layers"]
            / (k * a)
        ),
        consumer=s_injection(),
        slots=slots,
        consumer_Z_masks=[str(p) for p in consumer],
    )
