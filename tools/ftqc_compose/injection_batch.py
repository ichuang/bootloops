"""Exact accepted batch distributions and an ideal logical S consumer."""

from collections import defaultdict
from fractions import Fraction
import json

from ftqc_extractors.injection_circuit import paths, run
from ftqc_audit.injection_records import parse, responses


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


def coordinates(locations):
    if locations == "all":
        return None
    if not isinstance(locations, list):
        raise ValueError("locations must be all or a list of [operation, sites]")
    keys = set()
    for entry in locations:
        if (
            not isinstance(entry, list)
            or len(entry) != 2
            or type(entry[0]) is not int
            or not isinstance(entry[1], list)
            or any(type(q) is not int for q in entry[1])
        ):
            raise ValueError("invalid physical channel coordinate")
        key = (entry[0], tuple(entry[1]))
        if key in keys:
            raise ValueError("repeated channel")
        keys.add(key)
    return keys


def result(s):
    try:
        import stim
    except ImportError as exc:
        raise ValueError("injection-batch requires Stim") from exc

    if set(s) != {"kind", "case", "locations", "frame", "slots", "state_budget"}:
        raise ValueError("unsupported batch contract fields")
    case = s["case"]
    source, path, _, _ = paths(case)
    spec = json.loads(source.read_text())
    evidence = run(dict(case=case))
    c = stim.Circuit.from_file(path)
    parsed = parse(c)
    width = evidence["postselect_width"]
    k = c.num_observables
    frame = s["frame"]
    if (
        not isinstance(frame, list)
        or len(frame) != k
        or any(
            type(f) is not int or f < 0 or f.bit_length() > c.num_detectors
            for f in frame
        )
    ):
        raise ValueError("one detector parity mask per output frame bit")
    slots = s["slots"]
    if (
        not isinstance(slots, list)
        or not slots
        or any(type(i) is not int or not 0 <= i < k for i in slots)
        or len(set(slots)) != len(slots)
    ):
        raise ValueError("the consumer requires distinct output slots")
    budget = s["state_budget"]
    if type(budget) is not int or not 1 <= budget <= 2000000:
        raise ValueError("invalid finite distribution budget")
    selected = coordinates(s["locations"])
    groups = {}
    for op, sites, mode, p, v in responses(c, parsed):
        key = op, sites
        if selected is not None and key not in selected:
            continue
        d = v & ((1 << c.num_detectors) - 1)
        out = v >> c.num_detectors
        for i, f in enumerate(frame):
            out ^= ((d & f).bit_count() % 2) << i
        signature = (d & ((1 << width) - 1)) | (out << width)
        group = groups.setdefault(key, defaultdict(Fraction))
        group[signature] += p
    if selected is not None and groups.keys() != selected:
        raise ValueError("unknown physical channel")
    mass = {0: Fraction(1)}
    for group in groups.values():
        remainder = 1 - sum(group.values())
        if remainder < 0:
            raise ValueError("physical channel probabilities exceed one")
        group[0] += remainder
        if group[0] == 1:
            continue
        nxt = defaultdict(Fraction)
        for current, a in mass.items():
            for displacement, b in group.items():
                if b:
                    nxt[current ^ displacement] += a * b
                    if len(nxt) > budget:
                        raise ValueError(
                            "full joint distribution exceeds state budget; no guarantee issued"
                        )
        mass = nxt
    accepted = [mass.get(i << width, 0) for i in range(1 << k)]
    a = sum(accepted)
    if not a:
        raise ValueError("zero acceptance")
    distribution = [p / a for p in accepted]
    consumer = [Fraction(0)] * (1 << len(slots))
    for m, p in enumerate(distribution):
        index = sum(((m >> slot) & 1) << j for j, slot in enumerate(slots))
        consumer[index] += p
    return dict(
        status="proved",
        kind="accepted simultaneous Y source and ideal S consumers",
        scope="independent listed physical channels; all unlisted locations ideal; ideal terminal decoding and transfer",
        noise_scope="all physical channels"
        if selected is None
        else "selected physical channels",
        circuit_sha256=evidence["circuit_sha256"],
        source_sha256=evidence["source_sha256"],
        physical_channels=len(groups),
        raw_sites=spec["preparation"]["carriers"],
        frame=frame,
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
