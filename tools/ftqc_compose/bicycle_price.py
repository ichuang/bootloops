"""Exact full-attempt pricing with independent named physical record faults."""

from fractions import Fraction
from itertools import product
from random import Random
from ftqc_extractors.__main__ import digest
from ftqc_extractors.bicycle import mapping
from ftqc_faults.audit.injection_statistics import interval
from .bicycle_install import install
from .supports import fields


def convolve(left, right):
    result = {}
    for (mask, p), (other, q) in product(left.items(), right.items()):
        if p * q:
            result[mask ^ other] = result.get(mask ^ other, 0) + p * q
    return result


def probability(v):
    if type(v) not in (int, str) or not 0 <= Fraction(v) <= 1:
        raise ValueError("expected an exact probability")
    return Fraction(v)


def price(s):
    fields(
        s, "kind", "installation", "faults", *(["sampling"] if "sampling" in s else [])
    )
    study = s["installation"]
    if study["kind"] != "bb144-rm15":
        raise ValueError("expected an installed BB144 factory")
    p = install(study, True)
    m = mapping()
    words = [
        sum(
            ((x & r["physical_z"]).bit_count() % 2) << j
            for j, (x, z) in enumerate(m["basis"][:5])
        )
        for r in m["rotations"]
    ]
    signatures = {0: Fraction(1)}
    for row in study["source"]["inputs"]:
        src = row["receipt"]["result"]
        group = {}
        for key, value in src["accepted_mass"].items():
            delta = 0
            for member, column in row["ports"]:
                if int(key) >> member & 1:
                    delta ^= words[column]
            group[delta] = group.get(delta, 0) + Fraction(value)
        signatures = convolve(signatures, group)
    a = sum(signatures.values(), Fraction(0))
    if (
        a != Fraction(p["factory"]["source_accept"])
        or signatures.get(0, 0) + signatures.get(1, 0)
        != Fraction(p["factory"]["final_accept"])
        or signatures.get(1, 0) != Fraction(p["factory"]["accepted_wrong"])
    ):
        raise ValueError(
            "noiseless installation differs from native source and factory"
        )
    seen = set()
    witnesses = []
    channels = []
    steps = {step["label"]: step for step in p["steps"]}
    initial = dict(signatures)
    for fault in s["faults"]:
        idx = fault["column"]
        if type(idx) is not int or not 0 <= idx < 15 or idx in seen:
            raise ValueError("invalid or duplicate site")
        loc = p["fault_locations"][idx]
        if {k: v for k, v in fault.items() if k != "p"} != loc:
            raise ValueError("unknown physical location")
        seen.add(idx)
        prob = probability(fault["p"])
        op = steps[loc["step"]]
        template = p["templates"][op["arms"][0]]
        if template["ticks"][loc["tick"]] != [["M", loc["wire"]]] or template[
            "readouts"
        ][0]["records"] != [0]:
            raise ValueError("premeasurement X does not flip the named record")
        # This fault exchanges the two meter outcomes with no subsequent data
        # interaction. Compare the physical adaptive arms for both reports b.
        frame = steps[f"r{idx}.frame"]
        for branch in p["ideal"]["branches"]:
            key = sum(v << j for j, v in enumerate(branch["outcomes"]))
            if frame["arms"][key] == frame["arms"][key ^ 2]:
                raise ValueError("record flip did not change the physical Pauli frame")
        delta = words[idx]
        channels.append((delta, prob))
        witnesses.append(
            dict(
                location=loc,
                absolute_tick=op["start"] + loc["tick"],
                flipped_record=loc["record"],
                gross_detector_flips=[],
                final_check_mask=delta >> 1,
                output_Z=delta & 1,
            )
        )
        signatures = convolve(signatures, {0: 1 - prob, delta: prob})
    wrong = signatures.get(1, Fraction(0))
    final = wrong + signatures.get(0, Fraction(0))

    def ratio(x, y):
        return str(x / y) if y else None

    prep = p["preparation_end"]
    core = p["elapsed"] - prep
    fwidth = p["factory_width"]
    width = p["width"]
    result = dict(
        status="proved",
        schedule_sha256=digest(p["steps"]),
        source_accept=str(a),
        final_accept=str(final),
        accepted_wrong=str(wrong),
        conditional_output_Z=ratio(wrong, final),
        witnesses=witnesses,
        factory_only=dict(
            condition="source group available",
            peak_reserved_qubits=fwidth,
            attempt_ticks=core,
            expected_qubit_cycles_per_output=ratio(fwidth * core * a, final),
        ),
        inclusive=dict(
            peak_reserved_qubits=width,
            source_reject_ticks=prep,
            full_attempt_ticks=prep + core,
            expected_ticks_per_group=str(prep + a * core),
            expected_qubit_cycles_per_output=ratio(width * (prep + a * core), final),
        ),
        premises="independent repeated groups; early release on source rejection; full footprint reserved until release; named premeasurement X faults and carrier noise only; all other gates, idles and decoding ideal",
        whole_program_fault_distance="open",
    )
    if "sampling" in s:
        cfg = s["sampling"]
        fields(cfg, "shots", "seed")
        n = cfg["shots"]
        if type(n) is not int or not 1 <= n <= 10000000 or type(cfg["seed"]) is not int:
            raise ValueError("invalid sample settings")
        rng = Random(cfg["seed"])
        totals = sorted(initial.items())
        accepted = errors = 0
        for _ in range(n):
            u = rng.random()
            v = 0.0
            for mask, mass in totals:
                v += float(mass)
                if u < v:
                    break
            else:
                continue
            for delta, prob in channels:
                if rng.random() < float(prob):
                    mask ^= delta
            if mask < 2:
                accepted += 1
                errors += mask
        result["sample"] = dict(
            cfg,
            accepted=accepted,
            wrong=errors,
            acceptance_interval=interval(accepted, n),
            conditional_error_interval=interval(errors, accepted),
        )
    return result
