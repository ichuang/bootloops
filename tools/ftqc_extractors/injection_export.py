"""Replayed accepted source batches for distillation studies."""

from collections import defaultdict
from fractions import Fraction
import json

from .injection_circuit import paths, run
from ftqc_faults.audit.injection_records import parse, responses


from .algebra import Limit


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


def y_law(s):
    try:
        import stim
    except ImportError as exc:
        raise ValueError("injection-batch requires Stim") from exc

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
                        raise Limit("joint distribution exceeds state budget")
        mass = nxt
    accepted = [mass.get(i << width, 0) for i in range(1 << k)]
    a = sum(accepted)
    return spec, evidence, accepted, len(groups)


def named(value):
    if not isinstance(value, str) or not value.strip():
        raise ValueError("expected a nonempty identifier")
    return value


def export(s, budget):
    from .__main__ import digest

    attempt, block = named(s["attempt"]), named(s["block"])
    m = s["model"]
    evidence_type = m["type"]
    if evidence_type == "closed-Y":
        if set(m) != {"type", "case", "locations", "frame"}:
            raise ValueError("invalid closed Y model")
        spec, proof, accepted, count = y_law(dict(m, state_budget=budget))
        law = dict(
            evidence="exact-Y-sign",
            accepted_mass={str(i): str(v) for i, v in enumerate(accepted) if v},
            acceptance=str(sum(accepted)),
            locations=m["locations"],
            physical_channels=count,
            frame=m["frame"],
            acceptance_bits=proof["postselect_width"],
            circuit_sha256=proof["circuit_sha256"],
            preparation_sha256=proof["source_sha256"],
            target="Y",
        )
    elif evidence_type == "raw-T":
        from .injection_t import accepted_t
        spec, law = accepted_t(m,budget)
    elif evidence_type == "sampled-Y":
        from ftqc_faults.audit.injection_statistics import summarize

        if set(m) != {"type", "case", "sample"}:
            raise ValueError("invalid sampled source")
        proof = run(dict(case=m["case"]))
        spec = json.loads(paths(m["case"])[0].read_text())
        sample = m["sample"]
        if (
            sample["circuit_sha256"] != proof["circuit_sha256"]
            or sample["postselected_detectors"] != proof["postselect_width"]
        ):
            raise ValueError("sample has wrong circuit or acceptance")
        stats = summarize(
            sample,
            proof["outputs"],
            proof["cost"]["physical_qubits"],
            proof["cost"]["two_qubit_layers"],
        )
        law = dict(
            evidence="sampled-Y-sign",
            target="Y",
            sample=sample,
            statistics=stats,
            preparation_sha256=proof["source_sha256"],
            circuit_sha256=proof["circuit_sha256"],
        )
    else:
        raise ValueError("unknown injection evidence")
    return dict(
        schema="ftqc-injection-batch-1",
        status="proved",
        attempt=attempt,
        block=block,
        slots=list(range(len(spec["logicals"]))),
        signed_pairs=spec["logicals"],
        model_sha256=digest(m),
        spectator_scope="unspecified logicals of a freshly prepared block",
        **law,
    )
