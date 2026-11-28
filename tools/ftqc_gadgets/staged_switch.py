"""Exact CSS gauge fixing on the Steane block and its extended eight-site ancilla.

Word amplitudes are integers with one common normalization per branch. Both
logical basis columns must carry the same amplitude and branch probability.
"""

from ftqc_codes.codes import steane, rm15
from ftqc_codes.gf2 import dot, nullspace
from .staged import Counterexample, endpoint


def span(rows):
    words = [0]
    for row in rows:
        words += [word ^ row for word in words]
    return words


def syndrome(word, checks):
    return sum(dot(word, check) << i for i, check in enumerate(checks))


def geometry():
    small, large = steane(), rm15()
    low = small.hx
    high = [row << 8 for row in low] + [255 << 7]
    split = low + high
    # Dual gauge Zs pair individually with the three low Xs. Choose the
    # lightest matching RM stabilizer in each class, breaking ties by mask.
    zwords = sorted(span(large.hz), key=lambda v: (v.bit_count(), v))
    dual = [next(z for z in zwords if syndrome(z, low) == 1 << i) for i in range(3)]
    return small, large, split, nullspace(split + [127], 15), dual


def example(operation):
    small, _, _, _, dual = geometry()
    s = dict(kind="staged-gadget", operation=operation, blocks=["data"])
    if operation == "rm15-phase":
        return dict(s, exponents=[-1] * 15)
    if operation not in ("raise", "lower"):
        raise ValueError("unknown staged operation")
    up = operation == "raise"
    return dict(
        s, checks=dual if up else small.hx, corrections=span(small.hx if up else dual)
    )


def check(spec):
    op = spec["operation"]
    base = example(op)
    if set(spec) != set(base) or len(spec["blocks"]) != 1:
        raise ValueError("invalid single-block stage")
    small, large, split, split_z, dual = geometry()
    names = spec["blocks"]
    source = "steane" if op == "raise" else "rm15"
    target = "steane" if op == "lower" else "rm15"
    result = dict(
        status="proved",
        inputs=endpoint(source, names),
        outputs=endpoint(target, names),
        action="T" if op == "rm15-phase" else "I",
    )
    if op == "rm15-phase":
        powers = spec["exponents"]
        if len(powers) != 15 or any(
            type(p) is not int or p not in (-1, 1) for p in powers
        ):
            raise ValueError("expected fifteen signed native T gates")
        for bit in (0, 1):
            for w in span(large.hx):
                word = w ^ (32767 * bit)
                phase = sum(p for j, p in enumerate(powers) if word >> j & 1) % 8
                if phase != bit:
                    raise Counterexample(
                        "signed logical T phase", bit=bit, word=word, phase=phase
                    )
        result["ideal"] = dict(coset_words=32, modulus=8)
    else:
        checks, corrections = spec["checks"], spec["corrections"]
        if (
            len(checks) != 3
            or len(corrections) != 8
            or any(
                type(v) is not int or not 0 <= v < 32768 for v in checks + corrections
            )
        ):
            raise ValueError(
                "expected three check masks and eight record-indexed corrections"
            )
        # Reject substitutions before using the canonical gauge fault model.
        if checks != base["checks"]:
            raise ValueError("unsupported gauge generators")
        up = op == "raise"
        for bit in (0, 1):
            split_words = {w ^ (127 * bit) for w in span(split)}
            rm_words = {w ^ (32767 * bit) for w in span(large.hx)}
            for record, correction in enumerate(corrections):
                if up:
                    output = {
                        w ^ correction
                        for w in split_words
                        if syndrome(w, checks) == record
                    }
                    correct = output == rm_words
                else:
                    # Expand the three commuting X projectors and apply the Z frame.
                    amps = {}
                    for subset, x in enumerate(span(checks)):
                        for w in rm_words:
                            v = w ^ x
                            amps[v] = amps.get(v, 0) + (-1) ** (
                                dot(subset, record) ^ dot(v, correction)
                            )
                    correct = set(amps) == split_words and set(amps.values()) == {1}
                if not correct:
                    raise Counterexample(
                        "corrected branch isometry", record=record, logical_bit=bit
                    )
        result["ideal"] = dict(
            branches=8,
            branch_probability="1/8",
            ancillary_x=split[3:],
            split_z=split_z,
            discarded_ancilla="fixed pure CSS state",
        )
    from .staged_switch_faults import check_faults

    result["faults"] = check_faults(spec, (small, large, split, split_z, dual))
    return result
