"""Replay the noisy gauge measurements and derive an endpoint recovery."""

from collections import Counter
from ftqc_codes.gf2 import dot
from ftqc_decoders.staged import Recovery
from ftqc_faults.staged import cat_faults
from .staged import Counterexample


def check_faults(spec, geometry, rule):
    small, large, split, split_z, _ = geometry
    op = spec["operation"]
    up = op == "raise"
    hx, hz = (split, split_z) if op == "lower" else (large.hx, large.hz)
    decoder = Recovery(hx, hz, rule)
    counts = Counter()
    locations = set()

    def observe(history, x, z, location):
        witness = decoder.observe([], x, z, dict(location=location, branch=history))
        if witness:
            raise Counterexample("ambiguous endpoint recovery", **witness)

    observe([], 0, 0, "clean")
    if op == "rm15-phase":
        # An incoming one-site error conjugated by T lies in the span of the
        # four Paulis on that same site. The ideal syndrome distinguishes their
        # cosets, including the identity; checking the span covers coherence.
        for q in range(15):
            for x, z in ((1, 0), (0, 1), (1, 1)):
                observe([], x << q, z << q, ["T", q, x, z])
                counts["native-T"] += 1
                counts["incoming"] += 1
            locations.add(("native-T", q))
    else:
        axis = "Z" if up else "X"
        checks = spec["checks"]
        schedule = [g for g in range(3) for _ in range(3)]

        def replay(at, mode, x, z, flip):
            delta = x if up else z
            shifts = [
                flip if k == at else int(k > at) * dot(delta, checks[g])
                for k, g in enumerate(schedule)
            ]
            mask = sum((sum(shifts[3 * g : 3 * g + 3]) >= 2) << g for g in range(3))
            for record in range(8):
                correction = (
                    spec["corrections"][record] ^ spec["corrections"][record ^ mask]
                )
                observe(
                    record,
                    x ^ (correction if up else 0),
                    z ^ (0 if up else correction),
                    [at, mode, mask],
                )

        for q in range(15):
            kind = "ancilla-output" if up and q >= 7 else "incoming"
            for x, z in ((1, 0), (0, 1), (1, 1)):
                replay(-1, [kind, q, x, z], x << q, z << q, 0)
                counts[kind] += 1
            if kind != "incoming":
                locations.add((kind, q))
        for at, g in enumerate(schedule):
            for mode, x, z, flip in cat_faults(15, checks[g], axis):
                replay(at, mode, x, z, flip)
                counts[mode[0]] += 1
                locations.add((at, mode[0], mode[1]))
    return dict(
        order=1,
        alternatives=dict(counts),
        noisy_locations=len(locations),
        recovery=decoder.rows(),
        terminal="ideal full syndrome and recovery",
        ideal_resources=(
            []
            if op == "rm15-phase"
            else ["verified cat seeds", "classical processing and frame update"]
            + (
                ["extended eight-site CSS seed before local output faults"]
                if up
                else []
            )
        ),
    )
