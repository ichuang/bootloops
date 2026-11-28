"""Upper bounds for XOR syndromes under a local-stochastic support law."""

from fractions import Fraction

from .stochastic import probability
from .supports import fields


def _masks(values, bits):
    if not isinstance(values, list) or any(
        type(mask) is not int or mask < 0 or mask.bit_length() > bits
        for mask in values
    ):
        raise ValueError("expected a list of unsigned syndrome masks fitting bits")
    return set(values)


def syndrome_envelope(bits, classes, p):
    """Per-syndrome upper weights for disjoint physical locations. See GUIDE.md.

    Each class gives its location count and the possible XOR effects of one
    faulty location. Absent syndromes have weight zero. Weights can exceed one.
    """
    if type(bits) is not int or bits < 0 or not isinstance(classes, list):
        raise ValueError("expected a nonnegative syndrome width and a list of classes")
    rate = probability(p)
    scale = 10**48
    mass = {0: scale}
    for group in classes:
        fields(group, "count", "effects")
        count = group["count"]
        if type(count) is not int or count < 0:
            raise ValueError("location count must be a nonnegative integer")
        effects = _masks(group["effects"], bits) - {0}
        if not effects or not rate:
            continue
        for _ in range(count):
            # Use the previous table throughout: one site contributes at most one effect.
            added = {}
            for mask in effects:
                for source, weight in mass.items():
                    dest = source ^ mask
                    added[dest] = added.get(dest, 0) + weight
            for dest, weight in added.items():
                delta = -(-weight * rate.numerator // rate.denominator)
                mass[dest] = mass.get(dest, 0) + delta
    return {mask: Fraction(weight, scale) for mask, weight in mass.items()}


def syndrome_result(spec):
    fields(spec, "kind", "bits", "classes", "p", "event")
    bounds = syndrome_envelope(spec["bits"], spec["classes"], spec["p"])
    event = _masks(spec["event"], spec["bits"])
    upper = min(Fraction(1), sum((bounds.get(s, 0) for s in event), Fraction(0)))
    return {
        "status": "conditional",
        "syndrome_event_upper": str(upper),
        "premises": [
            "Each count denotes distinct physical locations, disjoint across classes.",
            "Every specified k-location subset is simultaneously faulty with probability at most p^k.",
            "The syndrome starts at zero and is the XOR of the selected effects; each faulty site may also contribute zero.",
        ],
    }
